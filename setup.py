import setuptools
import glob

with open('README.md', 'r') as fh:
    long_description = fh.read()

# Update version from VERSION file into module
with open('VERSION', 'r') as fversion:
    version = fversion.readline().rstrip()
with open('_version.py', 'wt') as fversion:
    fversion.write('__version__ = "'+version+'"')


name = 'RNAseqpipe'
setuptools.setup(
    name=name, # Replace with your own username
    version=version,
    author='You Duan',
    author_email='duanyou@outlook.com',
    description='Analysis pipeline for RNA-seq data',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/yduanBioinfo/RNAseqpipe',
    install_requires = ['dypylib>=0.0.3', 'pyper'],
    #packages=setuptools.find_packages(),
    package_dir={name: '.'},
    packages=[name] + ['.'.join((name, x)) for x in setuptools.find_packages()],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
    ],
    package_data = {
        name: ['confs/*'] + ['group_data/*'] + 
              glob.glob('test/**', recursive=True) + 
              glob.glob('**/*.R', recursive=True) + 
              glob.glob('**/*.sh', recursive=True)
    },
    scripts=[
        'expression/count_merge',
        'get_gene_length.py',
        'run_RNAseqpipe.py',
        'install/post_RNAseqpipe_install',
        'stat/merge_flagstat.py',
        'conv_clipping.sh'
    ],
    python_requires='>=3.6',
    include_package_data=True
)
