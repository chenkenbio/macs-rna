import struct
from MACS3.Signal.Prob import (poisson_cdf, factorial, chisq_logp_e,
                               binomial_cdf, binomial_cdf_inv, binomial_pdf)

def hexf(x):
    # bit-exact f64 representation
    return "%016x" % struct.unpack('Q', struct.pack('d', float(x)))[0]

rows = []
# poisson_cdf grid: (n, lam, lower, log10)
grid = []
ns = [0, 1, 5, 10, 50, 80, 100, 200, 500, 1000, 1500, 3000]
lams = [0.5, 1.0, 5.0, 50.0, 100.0, 500.0, 699.0, 700.0, 701.0, 1000.0, 2000.0]
for n in ns:
    for lam in lams:
        for lower in (0, 1):
            for log10 in (0, 1):
                try:
                    v = poisson_cdf(n, lam, bool(lower), bool(log10))
                except AssertionError:
                    continue
                grid.append(("poisson_cdf", n, "%r" % lam, lower, log10, hexf(v)))
for g in grid:
    rows.append("\t".join(str(x) for x in g))

# factorial
for n in [0, 1, 2, 5, 10, 20, 50, 100, 170]:
    rows.append("\t".join(["factorial", str(n), "0", "0", "0", hexf(factorial(n))]))

# chisq_logp_e: (x, df, log10)
chisq = [(10,2),(100,2),(1000,22),(10,4),(100,8),(1000,80),(54,6),(565,10),(7765,12)]
for (x, df) in chisq:
    for log10 in (0, 1):
        v = chisq_logp_e(float(x), df, bool(log10))
        rows.append("\t".join(["chisq_logp_e", str(x), str(df), "0", str(log10), hexf(v)]))

# binomial_cdf: (x, a, b, lower)
binc = [(20,1000,0.01),(200,1000,0.01),(0,100,0.5),(50,100,0.5),(100,100,0.5),(3,10,0.3)]
for (x,a,b) in binc:
    for lower in (0,1):
        v = binomial_cdf(x, a, b, bool(lower))
        rows.append("\t".join(["binomial_cdf", str(x), str(a), "%r"%b, str(lower), hexf(v)]))

# binomial_pdf: (x, a, b)
for (x,a,b) in binc:
    v = binomial_pdf(x, a, b)
    rows.append("\t".join(["binomial_pdf", str(x), str(a), "%r"%b, "0", hexf(v)]))

# binomial_cdf_inv: (cdf, a, b) -> int
for (cdf,a,b) in [(0.1,1000,0.01),(0.01,1000,0.01),(0.5,100,0.5),(0.99,100,0.5)]:
    v = binomial_cdf_inv(cdf, a, b)
    rows.append("\t".join(["binomial_cdf_inv", "%r"%cdf, str(a), "%r"%b, "0", str(int(v))]))

out = "/home/kenchen/Documents/macs-rna/rust/macs3-rs/crates/macs3-core/tests/fixtures/prob_golden.tsv"
with open(out,"w") as fh:
    fh.write("# func\targ1\targ2\targ3\targ4\tresult(f64 bits hex, except *_inv which is integer) -- generated with MACS3 3.0.4 python\n")
    fh.write("\n".join(rows)+"\n")
print("wrote", len(rows), "rows to", out)
