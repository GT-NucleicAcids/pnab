import os
from IPython.core.display import display, HTML


class view:

    def __init__(self, jsmol_path, molecule_path):
        jsmol_path = os.path.join(jsmol_path, '')

        self.html = '''<!doctype html>
<html>
<head>  
	<meta content="text/html; charset=UTF-8" http-equiv="content-type">
	<title>Minimal JSmol Demo</title>

	<!-- CSS Style Inline: -->
	<style type="text/css">
		/* Jmol Applet */
		/* defines height, width and orientation of the Jmol div element */
		#jmol_div{
			height: 500px;
			width:  500px;
			float: left;
		}		
	</style>

	<!-- Load Jmol javascript library -->
	<script type="text/javascript" src= "''' + jsmol_path + '''JSmol.min.js"></script>

	<!-- calls to jQuery and Jmol (inline) -->
	<script type="text/javascript">
		// Jmol readyFunction 
		// is called when Jmol is ready
		jmol_isReady = function(applet) {
			document.title = (applet._id + " - Jmol " + Jmol.___JmolVersion)
			Jmol._getElement(applet, "appletdiv").style.border="1px solid blue"
		}

		// initialize Jmol Applet
		var myJmol = "myJmol";
        var JmolPath = "''' + jsmol_path + '''"
		var Info = {
			width:   "100%",
			height:  "100%",
			color:   "#000000", //black
			use:     "HTML5",
			j2sPath: JmolPath + "j2s",// this needs to point to where the j2s directory is.
			jarPath: JmolPath + "java",// this needs to point to where the java directory is.
			jarFile: "JmolAppletSigned.jar",
			debug:   false,
			readyFunction: jmol_isReady,
			script:  'load "''' + molecule_path + '''" ;', // on-load Jmol script
	 		allowJavaScript: false,
			disableJ2SLoadMonitor: true,
		}

		// jQuery ready functions
		// is called when page has been completely loaded
		$(document).ready(function() {
			$("#jmol_div").html(Jmol.getAppletHtml(myJmol, Info))
		})
		var lastPrompt=0;
	</script>

</head>
<body>
<h1>''' + molecule_path + '''</h1>

<!-- DIV in which Jmol Applet will be placed -->
<div id='jmol_div'></div>

</body>
</html>'''


    def view(self):
       display(HTML(self.html))
