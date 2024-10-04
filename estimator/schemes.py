from sage.all import oo
from .nd import stddevf, Binary, CenteredBinomial, DiscreteGaussian, SparseTernary, UniformMod
from .lwe_parameters import LWEParameters
from .ntru_parameters import NTRUParameters
from .sis_parameters import SISParameters

# NIST PQC Round 3 Finalists

#
# Kyber
#
#
# https://pq-crystals.org/kyber/data/kyber-specification-round3-20210804.pdf
# Table 1, Page 11, we are ignoring the compression
#
# https://eprint.iacr.org/2020/1308.pdf
# Table 2, page 27, disagrees on Kyber 512

Kyber512 = LWEParameters(
    n=2 * 256,
    q=3329,
    Xs=CenteredBinomial(3),
    Xe=CenteredBinomial(3),
    m=2 * 256,
    tag="Kyber 512",
)

Kyber768 = LWEParameters(
    n=3 * 256,
    q=3329,
    Xs=CenteredBinomial(2),
    Xe=CenteredBinomial(2),
    m=3 * 256,
    tag="Kyber 768",
)

Kyber1024 = LWEParameters(
    n=4 * 256,
    q=3329,
    Xs=CenteredBinomial(2),
    Xe=CenteredBinomial(2),
    m=4 * 256,
    tag="Kyber 1024",
)

#
# Saber
#
#
# https://www.esat.kuleuven.be/cosic/pqcrypto/saber/files/saberspecround3.pdf
# Table 1, page 11
#
# https://eprint.iacr.org/2020/1308.pdf
# Table 2, page 27, agrees

LightSaber = LWEParameters(
    n=2 * 256,
    q=8192,
    Xs=CenteredBinomial(5),
    Xe=UniformMod(8),
    m=2 * 256,
    tag="LightSaber",
)

Saber = LWEParameters(
    n=3 * 256,
    q=8192,
    Xs=CenteredBinomial(4),
    Xe=UniformMod(8),
    m=3 * 256,
    tag="Saber",
)

FireSaber = LWEParameters(
    n=4 * 256,
    q=8192,
    Xs=CenteredBinomial(3),
    Xe=UniformMod(8),
    m=4 * 256,
    tag="FireSaber",
)

#
# NTRU
#
#

NTRUHPS2048509Enc = NTRUParameters(
    n=508,
    q=2048,
    Xe=SparseTernary(2048 / 16 - 1),
    Xs=UniformMod(3),
    m=508,
    tag="NTRUHPS2048509Enc",
)

NTRUHPS2048677Enc = NTRUParameters(
    n=676,
    q=2048,
    Xs=UniformMod(3),
    Xe=SparseTernary(2048 / 16 - 1),
    m=676,
    tag="NTRUHPS2048677Enc",
)

NTRUHPS4096821Enc = NTRUParameters(
    n=820,
    q=4096,
    Xs=UniformMod(3),
    Xe=SparseTernary(4096 / 16 - 1),
    m=820,
    tag="NTRUHPS4096821Enc",
)

NTRUHRSS701Enc = NTRUParameters(
    n=700,
    q=8192,
    Xs=UniformMod(3),
    Xe=UniformMod(3),
    m=700,
    tag="NTRUHRSS701",
)

#
# Dilithium
#
#
# https://pq-crystals.org/dilithium/data/dilithium-specification-round3-20210208.pdf
# Table 1, Page 8

Dilithium2_MSIS_WkUnf = SISParameters(
    n=256*4,
    q=8380417,
    length_bound=350209,
    m=256*9,
    norm=oo,
    tag="Dilithium2_MSIS_WkUnf"
)

Dilithium2_MSIS_StrUnf = SISParameters(
    n=256*4,
    q=8380417,
    length_bound=380929,
    m=256*9,
    norm=oo,
    tag="Dilithium2_MSIS_StrUnf"
)

Dilithium3_MSIS_WkUnf = SISParameters(
    n=256*6,
    q=8380417,
    length_bound=724481,
    m=256*6*2,
    norm=oo,
    tag="Dilithium3_MSIS_WkUnf"
)

Dilithium3_MSIS_StrUnf = SISParameters(
    n=256*6,
    q=8380417,
    length_bound=1048576,
    m=256*6*2,
    norm=oo,
    tag="Dilithium3_MSIS_StrUnf"
)

Dilithium5_MSIS_WkUnf = SISParameters(
    n=256*8,
    q=8380417,
    length_bound=769537,
    m=256*8*2,
    norm=oo,
    tag="Dilithium5_MSIS_WkUnf"
)

Dilithium5_MSIS_StrUnf = SISParameters(
    n=256*8,
    q=8380417,
    length_bound=1048576,
    m=256*8*2,
    norm=oo,
    tag="Dilithium5_MSIS_StrUnf"
)

NISTPQC_R3 = (
    Kyber512,
    Kyber768,
    Kyber1024,
    LightSaber,
    Saber,
    FireSaber,
    NTRUHPS2048509Enc,
    NTRUHPS2048677Enc,
    NTRUHPS4096821Enc,
    NTRUHRSS701Enc,
)

#
# Falcon
#
#
# https://falcon-sign.info/falcon.pdf
# Table 3.3 (P. 51)

Falcon512_Unf = SISParameters(
    n=512,
    q=12289,
    length_bound=5833.9072,
    m=1024,
    norm=2,
    tag="Falcon512_Unf"
)

Falcon512_SKR = NTRUParameters(
    n=512,
    q=12289,
    Xs=DiscreteGaussian(4.0532),
    Xe=DiscreteGaussian(4.0532),
    m=512,
    ntru_type='circulant',
    tag="Falcon512_SKR"
)

Falcon1024_Unf = SISParameters(
    n=1024,
    q=12289,
    length_bound=8382.4081,
    m=2048,
    norm=2,
    tag="Falcon1024_Unf"
)

Falcon1024_SKR = NTRUParameters(
    n=1024,
    q=12289,
    Xs=DiscreteGaussian(2.866),
    Xe=DiscreteGaussian(2.866),
    m=1024,
    ntru_type='circulant',
    tag="Falcon1024_SKR"
)

# FrodoKEM
# https://frodokem.org/files/FrodoKEM-specification-20210604.pdf#page=24

Frodo640 = LWEParameters(
    n=640,
    q=2**15,
    Xs=DiscreteGaussian(2.8),
    Xe=DiscreteGaussian(2.8),
    m=640 + 16,
    tag="Frodo640",
)

Frodo976 = LWEParameters(
    n=976,
    q=2**16,
    Xs=DiscreteGaussian(2.3),
    Xe=DiscreteGaussian(2.3),
    m=976 + 16,
    tag="Frodo976",
)

Frodo1344 = LWEParameters(
    n=1344,
    q=2**16,
    Xs=DiscreteGaussian(1.4),
    Xe=DiscreteGaussian(1.4),
    m=1344 + 16,
    tag="Frodo1344",
)

# HES v1.1

HESv111024128error = LWEParameters(
    n=1024,
    q=2**27,
    Xs=DiscreteGaussian(3.0),
    Xe=DiscreteGaussian(3.0),
    m=1024,
    tag="HESv11error",
)

HESv111024128ternary = LWEParameters(
    n=1024,
    q=2**27,
    Xs=UniformMod(3),
    Xe=DiscreteGaussian(3.0),
    m=1024,
    tag="HESv11ternary",
)

HESv11 = (HESv111024128error, HESv111024128ternary)


# FHE schemes

# TFHE
# https://tfhe.github.io/tfhe/security_and_params.html
# - Key-Switching key (LWE)
TFHE630 = LWEParameters(
    n=630,
    q=2**32,
    Xs=Binary,
    Xe=DiscreteGaussian(stddev=2 ** (-15) * 2**32),
    tag="TFHE630",
)
# - Bootstrapping key (Ring-LWE)
TFHE1024 = LWEParameters(
    n=1024,
    q=2**32,
    Xs=Binary,
    Xe=DiscreteGaussian(stddev=2 ** (-25) * 2**32),
    tag="TFHE1024",
)

# CONCRETE default parameter set for 128-bit security, chosen for
# faster boolean circuit evalutation than the TFHE_LIB parameters.
# With these parameters, the probability of error is upper-bounded by 2^-25.
# https://github.com/zama-ai/concrete/blob/4209e3366e8eb889e83720de3dc03b85778d3cec/concrete-boolean/src/parameters/mod.rs#L83
# - Key-Switching key (LWE)
Concrete_TFHE586 = LWEParameters(
    n=586,
    q=2**32,
    Xs=Binary,
    Xe=DiscreteGaussian(stddev=2 ** (-13.4) * 2**32),
    tag="Concrete_TFHE586",
)
# - Bootstrapping key (Ring-LWE)
Concrete_TFHE512 = LWEParameters(
    n=512,
    q=2**32,
    Xs=Binary,
    Xe=DiscreteGaussian(stddev=2 ** (-24.8) * 2**32),
    tag="Concrete_TFHE512",
)

# https://eprint.iacr.org/2018/421.pdf
# Table 3, page 55

TFHE16_500 = LWEParameters(
    n=500,
    q=2**32,
    Xs=Binary,
    Xe=DiscreteGaussian(stddev=2.43 * 10 ** (-5) * 2**32),
    tag="TFHE16_500",
)

TFHE16_1024 = LWEParameters(
    n=1024,
    q=2**32,
    Xs=Binary,
    Xe=DiscreteGaussian(stddev=3.73 * 10 ** (-9) * 2**32),
    tag="TFHE16_1024",
)

# https://eprint.iacr.org/2018/421.pdf
# Table 4, page 55
TFHE20_612 = LWEParameters(
    n=612,
    q=2**32,
    Xs=Binary,
    Xe=DiscreteGaussian(stddev=2 ** (-15) * 2**32),
    tag="TFHE20_612",
)

TFHE20_1024 = LWEParameters(
    n=1024,
    q=2**32,
    Xs=Binary,
    Xe=DiscreteGaussian(stddev=2 ** (-26) * 2**32),
    tag="TFHE20_1024",
)

# FHEW
# https://eprint.iacr.org/2014/816.pdf
# page 14

FHEW = LWEParameters(
    n=500,
    q=2**32,
    Xs=Binary,
    Xe=DiscreteGaussian(stddev=2 ** (-15) * 2**32),
    tag="FHEW",
)

# SEAL

# v2.0
# https://www.microsoft.com/en-us/research/wp-content/uploads/2016/09/sealmanual.pdf
# Table 3, page 19

SEAL20_1024 = LWEParameters(
    n=1024,
    q=2**48 - 2**20 + 1,
    Xs=UniformMod(3),
    Xe=DiscreteGaussian(stddev=3.19),
    tag="SEAL20_1024",
)

SEAL20_2048 = LWEParameters(
    n=2048,
    q=2**94 - 2**20 + 1,
    Xs=UniformMod(3),
    Xe=DiscreteGaussian(stddev=3.19),
    tag="SEAL20_2048",
)

SEAL20_4096 = LWEParameters(
    n=4096,
    q=2**190 - 2**30 + 1,
    Xs=UniformMod(3),
    Xe=DiscreteGaussian(stddev=3.19),
    tag="SEAL20_4096",
)

SEAL20_8192 = LWEParameters(
    n=8192,
    q=2**383 - 2**33 + 1,
    Xs=UniformMod(3),
    Xe=DiscreteGaussian(stddev=3.19),
    tag="SEAL20_8192",
)

SEAL20_16384 = LWEParameters(
    n=16384,
    q=2**767 - 2**56 + 1,
    Xs=UniformMod(3),
    Xe=DiscreteGaussian(stddev=3.19),
    tag="SEAL20_16384",
)

# v2.2
# https://www.microsoft.com/en-us/research/wp-content/uploads/2017/06/sealmanual_v2.2.pdf
# Table 3, page 20

SEAL22_2048 = LWEParameters(
    n=2048,
    q=2**60 - 2**14 + 1,
    Xs=UniformMod(3),
    Xe=DiscreteGaussian(stddev=3.19),
    tag="SEAL22_2048",
)

SEAL22_4096 = LWEParameters(
    n=4096,
    q=2**116 - 2**18 + 1,
    Xs=UniformMod(3),
    Xe=DiscreteGaussian(stddev=3.19),
    tag="SEAL22_4096",
)

SEAL22_8192 = LWEParameters(
    n=8192,
    q=2**226 - 2**26 + 1,
    Xs=UniformMod(3),
    Xe=DiscreteGaussian(stddev=3.19),
    tag="SEAL22_8192",
)

SEAL22_16384 = LWEParameters(
    n=16384,
    q=2**435 - 2**33 + 1,
    Xs=UniformMod(3),
    Xe=DiscreteGaussian(stddev=3.19),
    tag="SEAL22_16384",
)

SEAL22_32768 = LWEParameters(
    n=32768,
    q=2**889 - 2**54 - 2**53 - 2**52 + 1,
    Xs=UniformMod(3),
    Xe=DiscreteGaussian(stddev=3.19),
    tag="SEAL22_32768",
)

# The following are not parameters of actual schemes
# but useful for benchmarking

# HElib
# https://eprint.iacr.org/2017/047.pdf
# Table 1, page 6
# 80-bit security

HElib80_1024 = LWEParameters(
    n=1024,
    q=2**47,
    Xs=SparseTernary(32),
    Xe=DiscreteGaussian(stddev=3.2),
    tag="HElib80_1024",
)

HElib80_2048 = LWEParameters(
    n=2048,
    q=2**87,
    Xs=SparseTernary(32),
    Xe=DiscreteGaussian(stddev=3.2),
    tag="HElib80_2048",
)

HElib80_4096 = LWEParameters(
    n=4096,
    q=2**167,
    Xs=SparseTernary(32),
    Xe=DiscreteGaussian(stddev=3.2),
    tag="HElib80_4096",
)

# 120-bit security

HElib120_1024 = LWEParameters(
    n=1024,
    q=2**38,
    Xs=SparseTernary(32),
    Xe=DiscreteGaussian(stddev=3.2),
    tag="HElib80_1024",
)

HElib120_2048 = LWEParameters(
    n=2048,
    q=2**70,
    Xs=SparseTernary(32),
    Xe=DiscreteGaussian(stddev=3.2),
    tag="HElib80_2048",
)

HElib120_4096 = LWEParameters(
    n=4096,
    q=2**134,
    Xs=SparseTernary(32),
    Xe=DiscreteGaussian(stddev=3.2),
    tag="HElib80_4096",
)


# Test parameters from CHHS
# https://eprint.iacr.org/2019/1114.pdf
# Table 4, page 18

CHHS_1024_25 = LWEParameters(
    n=1024,
    q=2**25,
    Xs=SparseTernary(32),
    Xe=DiscreteGaussian(stddev=stddevf(8)),
    tag="CHHS_1024_25",
)

CHHS_2048_38 = LWEParameters(
    n=2048,
    q=2**38,
    Xs=SparseTernary(32),
    Xe=DiscreteGaussian(stddev=stddevf(8)),
    tag="CHHS_2048_38",
)

CHHS_2048_45 = LWEParameters(
    n=2048,
    q=2**45,
    Xs=SparseTernary(32),
    Xe=DiscreteGaussian(stddev=stddevf(8)),
    tag="CHHS_2048_45",
)

CHHS_4096_67 = LWEParameters(
    n=4096,
    q=2**67,
    Xs=SparseTernary(32),
    Xe=DiscreteGaussian(stddev=stddevf(8)),
    tag="CHHS_4096_67",
)

CHHS_4096_82 = LWEParameters(
    n=4096,
    q=2**82,
    Xs=SparseTernary(32),
    Xe=DiscreteGaussian(stddev=stddevf(8)),
    tag="CHHS_4096_82",
)
