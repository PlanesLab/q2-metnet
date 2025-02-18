from setuptools import setup, find_packages

# Setup copied from q2-emperor
setup(
    name="q2-metnet",
    version="2023.0.1",
    packages=find_packages(),
    author="Francesco Balzerani, Telmo Blasco, Luis Vitores Valcarcel, Francisco J. Planes",
    author_email="fbalzerani@tecnun.es, tblasco@tecnun.es, lvalcarcel@tecnun.es, fplanes@tecnun.es",
    description="Package to contextualize table of taxonomical frequency of samples to metabolic reconstruction (AGORA or AGREDA) and extract features based on score of active reactions or subsystems",
    license='',
    url="https://qiime2.org",
    entry_points={
        'qiime2.plugins':
        ['q2-metnet=q2_metnet.plugin_setup:plugin']
    },
    zip_safe=False,
    package_data={
        'q2_metnet': [
            'data/*/*',
            '_clustermap/assets/index.html',
            '_pca/assets/index.html',
            '_boxplot/assets/index.html'
        ]
    }
)
