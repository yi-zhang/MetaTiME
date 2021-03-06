import setuptools

#with open("README.md", "r", encoding="utf-8") as fh:
#    long_description = fh.read()

setuptools.setup(
    name='metatime',
    version='1.1.0',
    author='Yi Zhang',
    author_email='yiz@ds.dfci.harvard.edu',
    description='Testing installation of metatime',
    long_description='test metatime install',
    long_description_content_type="text/markdown",
	url='https://github.com/yi-zhang/MetaTiME.git',
    license='MIT',
    packages=['metatime'],
	package_data = {
		'metatime': ['pretrained/mec/*txt']
		},
    install_requires=['pandas','scanpy','anndata','matplotlib','adjustText', 'leidenalg'],
)
