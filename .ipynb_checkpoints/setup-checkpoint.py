from setuptools import setup

requirements = [
    'pandas',
    'numpy',
    'msprime'
]

setup(
    name='vcfsim',
    version='1.0.0-alpha',
    packages=['vcfsim'],
    entry_points={
        'console_scripts': [
            'vcfsim=vcfsim.__main__:main'
        ]
    },
    url='https://github.com/Pie115/VCFSimulator-SamukLab',
    license='MIT',
    author='Paimon Goulart',
    author_email='paimongoulart@gmail.com',
    description="vcfsim",
    install_requires=requirements,
    keywords='pixy',
    classifiers=[
        'Programming Language :: Python :: 3.9.13'
    ]
)