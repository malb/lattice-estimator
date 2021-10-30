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
    Xs=NoiseDistribution.CentredBinomial(3),
    Xe=NoiseDistribution.CentredBinomial(2),
    m=4 * 256,
    tag="Kyber 512",
)

Kyber768 = LWEParameters(
    n=3 * 256,
    q=3329,
    Xs=NoiseDistribution.CentredBinomial(2),
    Xe=NoiseDistribution.CentredBinomial(2),
    m=3 * 256,
    tag="Kyber 768",
)

Kyber1024 = LWEParameters(
    n=4 * 256,
    q=3329,
    Xs=NoiseDistribution.CentredBinomial(2),
    Xe=NoiseDistribution.CentredBinomial(2),
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
    Xs=NoiseDistribution.CentredBinomial(5),
    Xe=NoiseDistribution.UniformMod(7),
    m=2 * 256,
    tag="LightSaber",
)

Saber = LWEParameters(
    n=3 * 256,
    q=8192,
    Xs=NoiseDistribution.CentredBinomial(4),
    Xe=NoiseDistribution.UniformMod(7),
    m=3 * 256,
    tag="Saber",
)

FireSaber = LWEParameters(
    n=4 * 256,
    q=8192,
    Xs=NoiseDistribution.CentredBinomial(3),
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
