 <?php 
	$FileNameBaseArray = array(
		"Source" => "crux-2.0.Source.tar.gz",
		"Linux32" => "crux-2.0.Linux.i686.zip",
		"Linux64" => "crux-2.0.Linux.x86_64.zip",
		"OSX" => "crux-2.0.Darwin.x86_64.zip",
		"Windows32" => "crux-2.0.Windows.ix86-pc.zip",
	);
	$downloadType = $_POST['downloadtype'];
	$filename = $FileNameBaseArray[$downloadType];
  	$filepath = "crux-2.0/" . $filename;
	header("Location: http://cruxtoolkit.sourceforge.net/download/" . $filepath);
	exit;
 ?> 
