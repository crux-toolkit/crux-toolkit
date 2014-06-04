#!/bin/bash
scp *.html charlesegrant@web.sourceforge.net:/home/project-web/cruxtoolkit/htdocs
scp *.png charlesegrant@web.sourceforge.net:/home/project-web/cruxtoolkit/htdocs
scp *.gif charlesegrant@web.sourceforge.net:/home/project-web/cruxtoolkit/htdocs
scp crux.html charlesegrant@web.sourceforge.net:/home/project-web/cruxtoolkit/htdocs/index.html
scp crux.css charlesegrant@web.sourceforge.net:/home/project-web/cruxtoolkit/htdocs
scp -r data charlesegrant@web.sourceforge.net:/home/project-web/cruxtoolkit/htdocs
scp -r demos charlesegrant@web.sourceforge.net:/home/project-web/cruxtoolkit/htdocs
scp -r download charlesegrant@web.sourceforge.net:/home/project-web/cruxtoolkit/htdocs
