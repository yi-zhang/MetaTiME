import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='metatime',
    version='1.3.0',
    author='Yi Zhang',
    author_email='wingsyiz@gmail.com',
	description='Beta MetaTiME: annotate TME scRNA cell states',
    long_description  = long_description,
    long_description_content_type="text/markdown",
	url='https://github.com/yi-zhang/MetaTiME.git',
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
    ],
    packages=['metatime'],
	package_data = {
		'metatime': ['pretrained/mec/*tsv']
		},
    python_requires='>=3.6',
    install_requires=['pandas','scanpy','anndata','matplotlib','adjustText', 'leidenalg', 'harmonypy', 'scipy','seaborn'],
)
