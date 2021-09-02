import matplotlib.pyplot as plt
import numpy as np


def fftIndgen(n):
    """
    [summary]
    
    [extended_summary]
    
    Args:
        n ([type]): [description]
    """
    a = np.arange(0, n // 2)
    b = np.arange(1, n // 2 + 1)
    b = b[::-1]
    b = [-i for i in b]
    return a + b


def gaussian_random_field(Pk=lambda k: k ** (-3.0), size=100):
    """
    gaussian_random_field [summary]
    
    [extended_summary]
    
    Args:
        Pk ([type], optional): [description]. Defaults to lambdak:k**(-3.0).
        size (int, optional): [description]. Defaults to 100.
    
    Returns:
        [type]: [description]
    """

    def Pk2(kx, ky):
        if kx == 0 and ky == 0:
            return 0.0
        return np.sqrt(Pk(np.sqrt(kx ** 2 + ky ** 2)))

    noise = np.fft.fft2(np.random.normal(size=(size, size)))
    amplitude = np.zeros((size, size))
    for i, kx in enumerate(fftIndgen(size)):
        for j, ky in enumerate(fftIndgen(size)):
            amplitude[i, j] = Pk2(kx, ky)
    return np.fft.ifft2(noise * amplitude)


for alpha in [-4.0, -3.0, -2.0]:
    out = gaussian_random_field(Pk=lambda k: k ** alpha, size=256)
    plt.figure()
    plt.imshow(out.real, cmap="jet")

plt.show()
