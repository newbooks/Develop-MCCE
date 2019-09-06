# How to write documentation

This documentation is written in mkdocs, which uses markdown language and publishes on [github pages](https://gunnerlab.github.io/Develop-MCCE/).

## Installation of Mkdocs

```bash
#!/usr/bin/env bash
pip install --upgrade pip
pip install bs4
pip install unicode
pip install mkdocs
pip install mkdocs-material
pip install pymdown-extensions
pip install markdown-blockdiag
pip install markdown-include
```

This will install mkdocs, theme, and necessary plugins for documentation of this project.

## Writing and editing documentation
All documentation files reside under folder wiki/

  * ```mkdocs.yml``` : configuration file, also defines the menu
  * ```docs/users/``` : location of documentation files for users
  * ```docs/developers/``` : location of documentation files for developers  
  * ```site/``` : site web pages published by command ```mkdocs gh-deploy```

Under wiki/ directory, run 
```
mkdocs serve
```
will bring up a web service, and point browser to [http://localhost:8000](http://localhost:8000) to enable viewing 
and debugging the documentation.

 
## Deploy the site
To deploy site, you can

  * create a pull request so that the GunnerLab repo administrator can merge and deploy for you.

or
 
  * deploy directly by running ```gh-deply``` under directory wiki/. You need GunnerLab repo member permission to do 
  this.


## Markdown languange references