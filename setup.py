
from setuptools import setup, find_packages

setup(name='envoMatch',
      version='1.0',
      description='Search for peptide isotopic envelopes in ms1 spectra.',
      author='Aaron Maurais',
      url='https://github.com/ajmaurais/envoMatch',
      classifiers=['Development Status :: 4 - Beta',
        'Intended Audience :: SCIENCE/RESEARCH',
        'Topic :: Software Development :: Build Tools',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        ],
      packages=find_packages(),
      package_dir={'envoMatch':'envoMatch'},
      python_requires='>=3.6.*',
      install_requires=['pyteomics>=4.4', 'lxml>=4.4', 'tqdm', 'sortedcontainers>=2.3', 'pandas>=1.0', 'matplotlib>=3.3'],
      entry_points={'console_scripts': ['envoMatch=envoMatch:main']},
)


