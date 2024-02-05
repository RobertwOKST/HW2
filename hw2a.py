from math import sqrt, pi, exp

def GPDF(args):
    x, mu, sig = args
    fx = (1 / (sig * sqrt(2 * pi))) * exp(-0.5 * ((x - mu) / sig) ** 2)
    return fx

def Simpson(fx, mu, sig, lhl, rhl):
    n = 20
    h = (rhl - lhl) / n


    area = (h / 3) * (fx(lhl, mu, sig) + 4 * sum(fx(lhl + i * h, mu, sig) for i in range(1, n, 2))
                     + 2 * sum(fx(lhl + i * h, mu, sig) for i in range(2, n, 2))
                     + fx(rhl, mu, sig))

    return area

def Probability(PDF, args, c, GT=True):
    mu, sig = args
    lhl = mu - 5 * sig
    rhl = c if GT else mu + 2 * sig
    if GT:
        p = Simpson(lambda x, mu, sig: PDF((x, mu, sig)), mu, sig, lhl, c)
    else:
        p = Simpson(lambda x, mu, sig: PDF((x, mu, sig)), mu, sig, lhl, rhl)
    return p

def main():
    """
    Integrates the Gaussian probability density function between
    a left hand limit = (mean - 5*stDev) to a right hand limit = (c).
    :return: Nothing to return, just print results to screen.
    """

    #Region Testing User input
    mean = float(input("Population mean? "))
    stDev = float(input("Standard deviation? "))
    c = float(input("c value? "))
    GT = input("Probability greater than c?: ").lower() in ["y", "yes", "true"]

    #Region Testing GPDF
    fx = GPDF((0, mean, stDev))
    print("GPDF(0, {}, {}) = {:.5f}".format(mean, stDev, fx))

    #Region Testing Simpson
    p = Simpson(lambda x, mu, sig: GPDF((x, mu, sig)), mean, stDev, -5 * stDev, 0)  # should return 0.5
    print("Simpson(GPDF, {}, {}, -5 * {}, 0) = {:.5f}".format(mean, stDev, stDev, p))

    #Region Testing Probability
    result = Probability(GPDF, (mean, stDev), c, GT)
    print("Probability(GPDF, ({}, {}), {}, {}) = {:.5f}".format(mean, stDev, c, GT, result))

    print("P(x {} {} | N({}, {})) = {:.5f}".format(">" if GT else "<", c, mean, stDev, result))

if __name__ == "__main__":
    main()