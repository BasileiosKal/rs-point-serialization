use blstrs::{G1Affine, G2Affine};
use groupy::CurveAffine;
use num_bigint::BigUint;
use num_traits::identities::{Zero, One};
use mcore::bls12381::{ecp::ECP, ecp2::ECP2, big::BIG, fp2::FP2};

const BLS12381_MODULUS_HEX: &str = "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab";
const BLS12381_BASE_FIELD_SIZE: usize = 48;


/// Wrapper trait for different curve points.     
/// BASE_FIELD_SIZE: number of bytes of base field element representation 
/// 
/// TODO: use associated types for BASE_FIELD_SIZE (parametrized the trait instead
///       of using associated types as a work around to issue #60551)
trait Point<'a, const BASE_FIELD_SIZE: usize> {
    // hex encoding of the prime order of the subgroup
    const MODULUS_HEX: &'a str;

    // cast the BASE_FIELD_SIZE to a const associated type
    const BASE_FIELD_SIZE: usize = BASE_FIELD_SIZE;

    // return the x coordinate as a vec of byte-arrays
    fn x_to_coordinates(&self) -> Vec<[u8; BASE_FIELD_SIZE]>;

    // return the y coordinate as a vec of byte-arrays
    fn y_to_coordinates(&self) -> Vec<[u8; BASE_FIELD_SIZE]>;

    // get the y coordinate form the x coordinate
    fn get_y_from_x(x: &[u8], s_bit: bool) -> Vec<u8>;

    // test if point is infinity
    fn is_infinity(&self) -> bool;

    // return the compressed x coordinate.
    // x = (x_1, x_2, ..., x_m)
    // x_compressed = x_1 || x_2 || ... || x_m
    fn x_to_compressed(&self) -> Vec<u8> {
        let x = self.x_to_coordinates();
        let mut x_compressed: Vec<u8> = Vec::new();
        for x_i in x {
            x_compressed.extend_from_slice(&x_i);
        }
        x_compressed
    }

    // return the x coordinates from the compressed representation.
    // x_compressed = x_1 || x_2 || ... || x_m
    // x_coordinates = (x_1, x_2, ..., x_m)
    fn x_from_compressed(x_compressed: Vec<u8>) -> Vec<[u8; BASE_FIELD_SIZE]> {
        let res: Vec<[u8; BASE_FIELD_SIZE]> = x_compressed.rchunks(BASE_FIELD_SIZE).map(
            |el| {
                let c: [u8; BASE_FIELD_SIZE] = el.try_into().unwrap();
                c
            }
        ).collect();
        res.to_vec()
    }

    // deserialize the modulus
    fn modulus() -> BigUint {
        BigUint::from_bytes_be(&hex::decode(Self::MODULUS_HEX).unwrap())
    }

    fn point_y_sign(&self) -> bool {
        let y = self.y_to_coordinates();
        let modulus = Self::modulus();
    
        for y_i in y {
            // make y_i a bigInt
            let y_i_scalar = BigUint::from_bytes_be(&y_i);

            if !y_i_scalar.is_zero() {
                if y_i_scalar > modulus - y_i_scalar.clone() - BigUint::one() {
                    return true
                }
                return false
            }
        }
        return false
    }

    fn encode(&self) -> Vec<u8> {
        // number of bits remaining unused by the modulus
        let rb = Self::MODULUS_HEX.len() * 4 - Self::modulus().bits() as usize;

        // 2 bytes for the metadata
        let mut metadata = 0u16;
 
        let c_bit = 0x0004 as u16;
        let mut s_bit = 0x0000 as u16;
        if self.point_y_sign() {
            s_bit = 0x0001 as u16;
        };

        let mut i_bit = 0x0000 as u16;
        if self.is_infinity() {
            i_bit = 0x0002 as u16;
        }

        metadata |= c_bit;
        metadata |= i_bit;
        metadata |= s_bit;
        metadata <<= 8-rb;

        // transform the metadata from u16 to [u8; 2]
        let metadata_bytes: [u8; 2] = metadata.to_be_bytes();

        // if a coordinate encoding has 3 bits to spare, use them to "save" the
        // 3 metadata bits. In any other case, use an additional byte. Initialize
        // the extra byte to 0 and remove it later if not needed.
        let mut compressed: Vec<u8> = vec![0u8];

        let mut x = self.x_to_compressed();
        compressed.append(&mut x);

        compressed[0] = metadata_bytes[0];
        compressed[1] |= metadata_bytes[1];

        // strip the "least significant" byte if it is 0
        if compressed[0] == 0u8 {
            compressed.remove(0);
        }

        compressed
    }

    // returns the x, y coordinates of the point
    fn decode(serialized_input:  &Vec<u8>) -> (Vec<u8>, Vec<u8>) {
        let mut serialized = serialized_input.clone();
        
        // get the bits remaining from the modulus to encode the metadata
        let rb = Self::MODULUS_HEX.len() * 4 - Self::modulus().bits() as usize;

        // if length(serialized) == BASE_FIELD_SIZE (<=> rb >=3) then all the metadata bits are included
        let mut metadata_bytes = [0u8; 2];

        if rb < 3 {
            metadata_bytes[0] = serialized[0];
            metadata_bytes[1] = serialized[1];

            // the first byte is only used for metadata
            serialized.remove(0);
        } else {
            metadata_bytes[1] = serialized[0];
        }

        // mask away the non metadata-bits
        let mut mask = 0x07;
        mask <<= 8-rb;
        metadata_bytes[1] &= mask;
        let mut metadata = u16::from_be_bytes(metadata_bytes);

        // mask away the metadata bits from the x coordinate encoding
        serialized[0] &= !mask;

        // move the metadata bits to the least significant bit positions
        metadata >>= 8 - rb;
        
        let c_bit = metadata & 0x0004 == 0x0004;
        let i_bit = metadata & 0x0002 == 0x0002;
        let s_bit = metadata & 0x0001 == 0x0001;

        // check for point at infinity
        if i_bit {
            //TODO: probably the least efficient way to do this
            if serialized.clone().into_iter().any(|el| el != 0u8) {
                panic!("Expected the point at infinity (0). Got not 0 point")
            }

            // S_bit must not be set
            if s_bit { panic!("Got point at infinity, with its sign set") }

            // derive the size of the coordinates from the size of the input
            let size = if c_bit {serialized.len()} else {serialized.len() / 2};
            let mut x: Vec<u8> = Vec::with_capacity(size);
            let mut y: Vec<u8> = Vec::with_capacity(size);
            x.resize(size, 0u8);
            y.resize(size, 0u8);
            return (x, y)
        }

        // Check that C_bit is in fact set
        if !c_bit{
            // TODO: should decode from uncompressed
            panic!("Uncompressed point decoding is not supported yet")
        }

        // get the y coordinate from x, based on the s_bit
        let x_bytes = serialized;
        let y_bytes = Self::get_y_from_x(&x_bytes, s_bit);
        (x_bytes.to_vec(), y_bytes.to_vec())
    }
}


// Implement the encoding/decoding rules for BLS12-381, G1 subgroup
impl<'b> Point<'b, BLS12381_BASE_FIELD_SIZE> for G1Affine {
    const MODULUS_HEX: &'b str = BLS12381_MODULUS_HEX;

    fn x_to_coordinates(&self) -> Vec<[u8; Self::BASE_FIELD_SIZE]> {
        let mut coordinate = Vec::new();
        coordinate.push(
            self.x().to_bytes_be()
        );
        coordinate
    }

    fn y_to_coordinates(&self) -> Vec<[u8; Self::BASE_FIELD_SIZE]> {
        let mut coordinate = Vec::new();
        coordinate.push(
            self.y().to_bytes_be()
        );
        coordinate
    }

    fn get_y_from_x(x: &[u8], s_bit: bool) -> Vec<u8> {
        let _point = ECP::new_bigint(&BIG::frombytes(&x), s_bit as isize);
        let mut buf = [0u8; 48];
        _point.gety().tobytes(&mut buf);
        buf.to_vec()
    }

    fn is_infinity(&self) -> bool {
        self.is_zero()
    }
}


// Implement the encoding/decoding rules for BLS12-381, G2 subgroup
impl<'b> Point<'b, BLS12381_BASE_FIELD_SIZE> for G2Affine {
    const MODULUS_HEX: &'b str = BLS12381_MODULUS_HEX;

    fn x_to_coordinates(&self) -> Vec<[u8; Self::BASE_FIELD_SIZE]> {
        let mut coordinate = Vec::new();
        let x = self.x();
        coordinate.push(x.c1().to_bytes_be());
        coordinate.push(x.c0().to_bytes_be());
        coordinate
    }

    fn y_to_coordinates(&self) -> Vec<[u8; Self::BASE_FIELD_SIZE]> {
        let mut coordinate = Vec::new();
        let y = self.y();
        coordinate.push(y.c1().to_bytes_be());
        coordinate.push(y.c0().to_bytes_be());
        coordinate
    }

    fn get_y_from_x(x: &[u8], s_bit: bool) -> Vec<u8> {
        // get the x coordinate as an element of Fp2
        let x_coordinates = Self::x_from_compressed(x.to_vec());
        if x_coordinates.len() != 2 {
            panic!("Malformed x coordinate")
        }

        let x_fp2 = FP2::new_bigs(
            &BIG::frombytes(&x_coordinates[0]),
            &BIG::frombytes(&x_coordinates[1])
        );

        // calculate the y coordinate given the s metadata bit
        let _point = ECP2::new_fp2(&x_fp2, !s_bit as isize); // TODO: Not sure why the negation is needed here. For some reason the opposite hint is used.
        let mut buf = [0u8; 2 * BLS12381_BASE_FIELD_SIZE];
        _point.gety().tobytes(&mut buf);

        buf.to_vec()
    }

    fn is_infinity(&self) -> bool {
        self.is_zero()
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use rand_core::SeedableRng;
    use rand_xorshift::XorShiftRng;
    use groupy::CurveProjective;
    use blstrs::{G1Projective, G2Projective, Fp2};

    fn _concat_bls12381_coordinates(s: Fp2) -> Vec<u8> {
        [
            s.c1().to_bytes_be(),
            s.c0().to_bytes_be()
        ].concat()
    }

    #[test]
    fn bls12381_g1_test() {
        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);
    
        let rand_point_affine: G1Affine = G1Projective::random(&mut rng).into_affine();
        
        // compressed point
        let mut encoded_point = rand_point_affine.encode();

        // native compression from blstrs
        let native_encoded_point = rand_point_affine.to_compressed().to_vec();
        assert_eq!(encoded_point, native_encoded_point);

        // un-compress point
        let point_coordinates = G1Affine::decode(&mut encoded_point);

        assert_eq!(point_coordinates.0, rand_point_affine.x().to_bytes_be().to_vec());
        assert_eq!(point_coordinates.1, rand_point_affine.y().to_bytes_be().to_vec());
    }

    #[test]
    fn bls12381_g2_test() {
        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0x62,
        ]);
    
        let mut rand_point = G2Projective::random(&mut rng);
        let rand_point_affine: G2Affine = rand_point.into_affine();
        
        // compressed point
        let mut encoded_point = rand_point_affine.encode();

        // native compression from blstrs
        let native_encoded_point = rand_point_affine.to_compressed().to_vec();
        assert_eq!(encoded_point, native_encoded_point);

        // test -P
        rand_point.negate();
        let neg_rand_point_affine = rand_point.into_affine();

        // compressed point
        let mut neg_encoded_point = neg_rand_point_affine.encode();

        // native compression from blstrs
        let neg_native_encoded_point = neg_rand_point_affine.to_compressed().to_vec();
        assert_eq!(neg_encoded_point, neg_native_encoded_point);

        // un-compress point
        let point_decoded = G2Affine::decode(&mut encoded_point);
        assert_eq!(point_decoded.0, _concat_bls12381_coordinates(rand_point_affine.x()));
        assert_eq!(point_decoded.1, _concat_bls12381_coordinates(rand_point_affine.y()));

        // un-compress negative point
        let neg_point_decoded = G2Affine::decode(&mut neg_encoded_point);
        assert_eq!(neg_point_decoded.0, _concat_bls12381_coordinates(neg_rand_point_affine.x()));
        assert_eq!(neg_point_decoded.1, _concat_bls12381_coordinates(neg_rand_point_affine.y()));
    }

    #[test]
    fn bls12381_infinity_encoding_test() {
        // G1
        let g1_point_zero: G1Affine = G1Projective::zero().into_affine();
        let encoded_g1_zero = g1_point_zero.encode();
        assert_eq!(encoded_g1_zero, g1_point_zero.to_compressed());

        let decoded_zero = G1Affine::decode(&encoded_g1_zero);
        assert_eq!(decoded_zero.0, g1_point_zero.x().to_bytes_be());
        assert_eq!(decoded_zero.1, g1_point_zero.y().to_bytes_be());

        // G2
        let g2_point_zero: G2Affine = G2Projective::zero().into_affine();
        let encoded_g2_zero = g2_point_zero.encode();
        assert_eq!(encoded_g2_zero, g2_point_zero.to_compressed());

        let decoded_g2_zero = G2Affine::decode(&encoded_g2_zero);
        assert_eq!(decoded_g2_zero.0, _concat_bls12381_coordinates(g2_point_zero.x()));
        assert_eq!(decoded_g2_zero.1, _concat_bls12381_coordinates(g2_point_zero.y()));
    }

}
