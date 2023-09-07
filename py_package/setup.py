from setuptools import setup, find_packages
 
classifiers = [
   'Development Status :: 5 - Production/Stable',
  'Operating System :: Microsoft :: Windows :: Windows 10',
  'Operating System :: POSIX :: Linux',
  'License :: OSI Approved :: MIT License',
  'Programming Language :: Python :: 3'
]
 
setup(
  name='SpatialPrompt',
  version='0.0.5',   #mandatory change
  description='Scalable and efficient cell type deconvolution and clustering in spatial transcriptomics',
  long_description=open('README.md').read() ,
  url='',  
  author='Asish Kumar Swain',
  author_email='swainkasish@gmail.com',
  license='MIT', 
  classifiers=classifiers,
  keywords='willupdate', 
  packages=find_packages(),
  install_requires=['pandas>=1.3.5',
        'numpy>=1.21.6',
        'scikit-learn>=1.0.2',
        "alive_progress>=3.0.1",
        "scipy>=1.7.3"] 
)
