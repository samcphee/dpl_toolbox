import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='dpl_toolbox',
    version='0.0.1',
    author='Scott McPhee',
    author_email='sahmcphee@gmail.com',
    description='DPL Toolbox',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/samcphee/dpl_toolbox',
    project_urls = {
        "Bug Tracker": "https://github.com/samcphee/dpl_toolbox/issues"
    },
    license='MIT',
    packages=['dpl_toolbox'],
    install_requires=['pandas', 'numpy', 'pyopenms','openpyxl'],
)
