from setuptools import setup, find_packages


with open("README.md", "r") as file: 
    long_description = file.read()

setup(
    name='ARBOLpy',
    version='0.0.8',
    packages=find_packages(),
    install_requires=[
        'scanpy>=1.9.0',
        'anytree>=2.8.0',
        'pandas>=1.0.5',
        'numpy>=1.18.5',
        'matplotlib>=3.5.1',
        'scipy>=1.4.1',
        'leidenalg>=0.8'
    ],
    python_requires=">=3.8",
    author='Kyle Kimler, Ruben van Esch',
    author_email='kkimler@broadinstitute.org',
    description='python implementation of ARBOL scRNAseq iterative tiered clustering\nhttps://github.com/jo-m-lab/ARBOL',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/jo-m-lab/ARBOLpy',
    project_urls={
        "Bug Tracker": 'https://github.com/jo-m-lab/ARBOLpy/issues', 
        "Source code": 'https://github.com/jo-m-lab/ARBOLpy'
    }
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
)
