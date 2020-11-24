from distutils.core import setup
setup(
        name='orsaytrace',
        packages = ['orsaytrace'],
        version = '1.0',
        license = 'MIT',
        description = 'Flexible ray tracing for optics',
        author = 'Yves Auad',
        author_email = 'yves.maia-auad@universite-paris-saclay.fr',
        url = 'https://github.com/yvesauad/OrsayTrace',
        download_url = 'https://github.com/yvesauad/OrsayTrace.git',
        keywords = ['optics', 'ray', 'tracing']
        install_requires = [
            'validators',
            'beautifulsoup4',
            ],
        classifiers = [
            'Development Status :: 3 - Alpha',
            'Intended Audience :: Developers',
            'Topic :: Software Development :: Build Tools',
            'License :: OSI Approved :: MIT License',
            'Programming Language :: Python :: 3.6',
            ],
        )

            
