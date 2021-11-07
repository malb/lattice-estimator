from .nd import NoiseDistribution
from .lwe import LWEParameters

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
    Xs=NoiseDistribution.CenteredBinomial(3),
    Xe=NoiseDistribution.CenteredBinomial(2),
    m=4 * 256,
    tag="Kyber 512",
)

Kyber768 = LWEParameters(
    n=3 * 256,
    q=3329,
    Xs=NoiseDistribution.CenteredBinomial(2),
    Xe=NoiseDistribution.CenteredBinomial(2),
    m=3 * 256,
    tag="Kyber 768",
)

Kyber1024 = LWEParameters(
    n=4 * 256,
    q=3329,
    Xs=NoiseDistribution.CenteredBinomial(2),
    Xe=NoiseDistribution.CenteredBinomial(2),
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
    Xs=NoiseDistribution.CenteredBinomial(5),
    Xe=NoiseDistribution.UniformMod(7),
    m=2 * 256,
    tag="LightSaber",
)

Saber = LWEParameters(
    n=3 * 256,
    q=8192,
    Xs=NoiseDistribution.CenteredBinomial(4),
    Xe=NoiseDistribution.UniformMod(7),
    m=3 * 256,
    tag="Saber",
)

FireSaber = LWEParameters(
    n=4 * 256,
    q=8192,
    Xs=NoiseDistribution.CenteredBinomial(3),
    Xe=NoiseDistribution.UniformMod(7),
    m=4 * 256,
    tag="FireSaber",
)

NTRUHPS2048509Enc = LWEParameters(
    n=508,
    q=2048,
    Xe=NoiseDistribution.SparseTernary(508, 2048 / 16 - 1),
    Xs=NoiseDistribution.UniformMod(3),
    m=508,
    tag="NTRUHPS2048509Enc",
)

NTRUHPS2048677Enc = LWEParameters(
    n=676,
    q=2048,
    Xs=NoiseDistribution.UniformMod(3),
    Xe=NoiseDistribution.SparseTernary(676, 2048 / 16 - 1),
    m=676,
    tag="NTRUHPS2048677Enc",
)

NTRUHPS4096821Enc = LWEParameters(
    n=820,
    q=4096,
    Xs=NoiseDistribution.UniformMod(3),
    Xe=NoiseDistribution.SparseTernary(820, 4096 / 16 - 1),
    m=820,
    tag="NTRUHPS4096821Enc",
)

NTRUHRSS701Enc = LWEParameters(
    n=700,
    q=8192,
    Xs=NoiseDistribution.UniformMod(3),
    Xe=NoiseDistribution.UniformMod(3),
    m=700,
    tag="NTRUHRSS701",
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

HESv111024128error = LWEParameters(
    n=1024,
    q=2 ** 27,
    Xs=NoiseDistribution.DiscreteGaussian(3.0),
    Xe=NoiseDistribution.DiscreteGaussian(3.0),
    m=1024,
    tag="HESv11error",
)

HESv111024128ternary = LWEParameters(
    n=1024,
    q=2 ** 27,
    Xs=NoiseDistribution.UniformMod(3),
    Xe=NoiseDistribution.DiscreteGaussian(3.0),
    m=1024,
    tag="HESv11ternary",
)

HESv11 = (HESv111024128error, HESv111024128ternary)
