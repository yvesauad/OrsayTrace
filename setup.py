#from distutils.core import setup
from setuptools import setup

proj_urls = {
                'Documentation': 'https://orsaytrace.readthedocs.io/en/latest/index.html',
                'Source Code': 'https://github.com/yvesauad/orsaytrace'
            }

setup(
        name='orsaytrace',
        packages = ['orsaytrace'],
        version = '1.1.11',
        license = 'MIT',
        description = 'Flexible ray tracing for optics',
        long_description = open('README.rst').read(),
        #long_description_content_type = 'text/markdown',
        author = 'Yves Auad',
        author_email = 'yves.maia-auad@universite-paris-saclay.fr',
        url = 'https://github.com/yvesauad/OrsayTrace',
        project_urls = proj_urls,
        download_url = 'https://github.com/yvesauad/OrsayTrace/archive/v1.1.3.tar.gz',
        keywords = ['optics', 'ray', 'tracing'],
        install_requires = [
            'numpy',
            'tqdm',
            'matplotlib',
            ],
        classifiers = [
            'Development Status :: 3 - Alpha',
            'Intended Audience :: Developers',
            'Topic :: Software Development :: Build Tools',
            'License :: OSI Approved :: MIT License',
            'Programming Language :: Python :: 3.9',
            ],
        )
