 <?php 
  # Set names of files available for download
   if (isset($_POST['release_button'])) {
     $FileNameBaseArray = array(
       "Source" => "crux-2.0.Source.tar.gz",
       "Linux32" => "crux-2.0.Linux.i686.zip",
       "Linux64" => "crux-2.0.Linux.x86_64.zip",
       "OSX" => "crux-2.0.Darwin.x86_64.zip",
       "Windows32" => "crux-2.0.Windows.ix86-pc.zip",
     );
     $directory = "crux-2.0/";
   }
    else {
      # File names for daily build include SVN revsion number.
      # Read revision number from file.
      $version = chop(file_get_contents("daily/latest-build.txt"));
      $FileNameBaseArray = array(
        "Source" => "crux-2.0.$version.Source.tar.gz",
        "Linux32" => "crux-2.0.$version.Linux.i686.zip",
        "Linux64" => "crux-2.0.$version.Linux.x86_64.zip",
        "OSX" => "crux-2.0.$version.Darwin.x86_64.zip",
        "Windows32" => "crux-2.0.$version.Windows.i686.zip",
      );
      $directory = "daily/";
    }
    # Redriect to file to be downloaded.
	  $downloadType = $_POST['downloadtype'];
	  $filename = $FileNameBaseArray[$downloadType];
    $filepath = $directory . $filename;
    header("Location: http://cruxtoolkit.sourceforge.net/download/" . $filepath);
    exit;
 ?> 
