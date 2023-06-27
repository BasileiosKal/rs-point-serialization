# EC Point Compressed Serialization PoC

PoC implementation of the elliptic curve point serialization procedure defined [here](https://github.com/BasileiosKal/ec-point-encoding)

## Build

This lib uses the rust crate of the [MIRACL library](https://github.com/miracl/core). Follow the steps [here](https://github.com/miracl/core/tree/master/rust#readme) to build the crate. Then add the path of the "mcore" file in the [cargo.toml](https://github.com/BasileiosKal/rs-point-serialization/blob/main/Cargo.toml#L10) file.