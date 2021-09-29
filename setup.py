from distutils.core import setup
setup(
  name = 'ou_Axion_limit',         # How you named your package folder (MyLib)
  packages = ['ou_Axion_limit'],   # Chose the same as "name"
  version = '0.5',      # Start with a small number and increase it with every change you make
  license='MIT',        # Chose a license from here: https://help.github.com/articles/licensing-a-repository
  description = 'Calculate the AXion G_gamma limit',   # Give a short description about your library
  long_description = open("CHANGELOG.txt").read(),
  author = 'OUYANG MINWEI',                   # Type in your name
  author_email = 'wesley91345@gmail.com',      # Type in your E-Mail
  url = 'https://github.com/OuYangMinOa/ou_Axion_limit',   # Provide either the link to your github or to your website
  download_url = '',    # I explain this later on
  keywords = ['Axion', 'TASEH', 'drak matter'],   # Keywords that define your package best
  install_requires=[            # I get to this in a second
          'matplotlib',
          'numpy',
      ],
  classifiers=[
    'Development Status :: 5 - Production/Stable',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
    'Intended Audience :: Developers',      # Define that your audience are developers
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',   # Again, pick a license
    'Programming Language :: Python :: 3',      #Specify which pyhton versions that you want to support
  ],
)