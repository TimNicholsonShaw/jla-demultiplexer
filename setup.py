from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

long_description = (here / 'README.md').read_text(encoding='utf-8')

setup(
    name='jla-demultiplexer', 
    version='0.0.5',
    description='Internal lab use tool for 3 end analysis',
    long_description=long_description,
    long_description_content_type='text/markdown', 
    url='https://github.com/TimNicholsonShaw/jla-demultiplexer',  
    author='Tim Nicholson-Shaw',
    author_email='timnicholsonshaw@gmail.com',
    classifiers=[ 
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3 :: Only',
    ],
    keywords='Sequencing, ncRNAs, Bioinformatics',
    package_dir={'': 'src'},
    packages=find_packages(where='src'),
    python_requires='>=3.6',
    install_requires=["jla-tailer>=0.1.5", "Affirmations>=0.0.10"],

    entry_points={
        'console_scripts': [
            'fastqbreakdown=deMultiplexer:fastqBreakDown',
            'demultiplexer=deMultiplexer:manifestAlign',
            'jla-trim=deMultiplexer:trimDeDup'
        ],
    },

    project_urls={ 
        'Lab Website': 'https://labs.biology.ucsd.edu/lykkeandersen/index.html',
    },
)