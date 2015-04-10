 <?php 
  # Set names of files available for download
   if (isset($_POST['release_button'])) {
     $FileNameBaseArray = array(
       "Source" => "crux-2.1.Source.tar.gz",
       "Linux32" => "crux-2.1.Linux.i686.zip",
       "Linux64" => "crux-2.1.Linux.x86_64.zip",
       "OSX" => "crux-2.1.Darwin.x86_64.zip",
       "Windows" => "crux-2.1.Windows.ix86-pc.zip",
     );
     $directory = "crux-2.1/";
   }
    else {
      # File names for daily build include SVN revsion number.
      # Read revision number from file.
      $version = chop(file_get_contents("daily/latest-build.txt"));
      $FileNameBaseArray = array(
        "Source" => "crux-2.1.$version.Source.tar.gz",
        "Linux32" => "crux-2.1.$version.Linux.i686.zip",
        "Linux64" => "crux-2.1.$version.Linux.x86_64.zip",
        "OSX" => "crux-2.1.$version.Darwin.x86_64.zip",
        "Windows" => "crux-2.1.$version.Windows.AMD64.zip",
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
