function logGamma(z) {

    const c = [
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313
    ];

    if (z < 0.5) {
        return Math.log(Math.PI)
             - Math.log(Math.sin(Math.PI * z))
             - logGamma(1 - z);
    }

    z -= 1;

    let x = 0.99999999999980993;

    return x;
}
