from distutils.core import setup
setup(
  name         = 'ou_Axion_limit',         # your package folder (ou_Axion_limit)
  packages     = ['ou_Axion_limit'],   # Chose the same as "name"
  version      = '1.0.0',      # your version
  license      = 'MIT',        # https://help.github.com/articles/licensing-a-repository
  description  = 'Calculate the AXion G_gamma limit',   # Give a short description about your library
  long_description = open("CHANGELOG.txt").read(),     # Give a long description about your library
  author       = 'OUYANG MINWEI',                    # your name
  author_email = 'wesley91345@gmail.com',      # your E-Mail
  url          = 'https://github.com/OuYangMinOa/ou_Axion_limit',   #  your github or to your website
  download_url = '',
  keywords     = ['Axion', 'TASEH', 'drak matter'],   # Keywords that define your package best
  install_requires = [            # required package
          'matplotlib',
          'numpy',
      ],
  classifiers = [
    'Development Status :: 5 - Production/Stable',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
    'Intended Audience :: Developers',       # Define that your audience are developers
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',   # Again, pick a license
    'Programming Language :: Python :: 3',      #Specify which pyhton versions that you want to support
  ],
)