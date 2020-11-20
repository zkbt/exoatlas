from exoatlas import *

def test_reflection(telescope_name='JWST', wavelength=1*u.micron):
    with mock.patch('builtins.input', return_value=""):
        t = TransitingExoplanets()

    w = 1*u.micron
    fi, ax = plt.subplots(1, 2, figsize=(8, 4))
    for i, per_transit in enumerate([True, False]):
        BubblePanel(xaxis=StellarBrightnessTelescope(wavelength=wavelength, telescope_name=telescope_name),
                    yaxis=Reflection(wavelength=w),
                    size=ReflectionSNR(wavelength=wavelength, telescope_name=telescope_name, per_transit=per_transit),
                    size_normalization=10).build(t, ax=ax[i])
        d = {True:'Transit', False:'Hour'}
        plt.title(f'S/N in 1 {d[per_transit]}')
    plt.tight_layout()


def test_emission(telescope_name='JWST', wavelength=5*u.micron):
    with mock.patch('builtins.input', return_value=""):
        t = TransitingExoplanets()

    w = 5*u.micron
    fi, ax = plt.subplots(1, 2, figsize=(8, 4))
    for i, per_transit in enumerate([True, False]):
        BubblePanel(xaxis=StellarBrightnessTelescope(wavelength=wavelength, telescope_name=telescope_name),
                    yaxis=Emission(wavelength=wavelength),
                    size=EmissionSNR(wavelength=wavelength, telescope_name=telescope_name),
                    size_normalization=10).build(t, ax=ax[i])
        d = {True:'Transit', False:'Hour'}
        plt.title(f'S/N in 1 {d[per_transit]}')
    plt.tight_layout()

def test_transmission(telescope_name='JWST', wavelength=5*u.micron):
    with mock.patch('builtins.input', return_value=""):
        t = TransitingExoplanets()

    w = 5*u.micron
    fi, ax = plt.subplots(1, 2, figsize=(8, 4))
    for i, per_transit in enumerate([True, False]):
        BubblePanel(xaxis=StellarBrightnessTelescope(wavelength=wavelength, telescope_name=telescope_name),
                    yaxis=Transmission(),
                    size=TransmissionSNR(wavelength=wavelength, telescope_name=telescope_name),
                    size_normalization=10).build(t, ax=ax[i])
        d = {True:'Transit', False:'Hour'}
        plt.title(f'S/N in 1 {d[per_transit]}')
    plt.tight_layout()

if __name__ == '__main__':
    outputs = {k.split('_')[-1]:v()
               for k, v in locals().items()
               if 'test_' in k}
