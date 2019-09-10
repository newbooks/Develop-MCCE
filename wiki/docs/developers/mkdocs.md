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

## Writing and editing
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

 
## Deploying
To deploy site, you can

  * create a pull request so that the GunnerLab repo administrator can merge and deploy for you.

or
 
  * deploy directly by running ```gh-deploy``` under directory wiki/. You need GunnerLab repo member permission to do 
  this.

If ```gh-deploy``` reports conflict in special branch gh-pages, use "--force" switch to overwrite:
```
mkdocs gh-deploy --force
```

The published site is on: [https://gunnerlab.github.io/Develop-MCCE/](https://gunnerlab.github.io/Develop-MCCE/).


## Markdown languange references

**Basic syntax:** [https://www.markdownguide.org/cheat-sheet/](https://www.markdownguide.org/cheat-sheet/)

**More on syntax:** [https://alinex.gitlab.io/env/mkdocs/](https://alinex.gitlab.io/env/mkdocs/)

**Diagram:** [http://blockdiag.com/en/blockdiag/](http://blockdiag.com/en/blockdiag/)