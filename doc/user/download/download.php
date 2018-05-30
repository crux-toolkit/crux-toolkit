<?php 
  function getRemoteIPAddress() {
    if (!empty($_SERVER['HTTP_CLIENT_IP'])) {
        return $_SERVER['HTTP_CLIENT_IP'];
    } else if (!empty($_SERVER['HTTP_X_FORWARDED_FOR'])) { 
        return $_SERVER['HTTP_X_FORWARDED_FOR'];
    }
    return $_SERVER['REMOTE_ADDR'];
  }

  function log_download($download_filename) {

     $log_filename = "downloads.txt";
   
    // Get time of request
    if( ($time = $_SERVER['REQUEST_TIME']) == '') {
      $time = time();
    }
   
    // Get IP address
    $remote_addr = getRemoteIPAddress();
    if($remote_addr  == '') {
      $remote_addr = "REMOTE_ADDR_UNKNOWN";
    }
   
    // Format the date and time
    $date = date("Y-m-d H:i:s", $time);
   
    // Append to the log file
    if($fd = @fopen($log_filename, "a")) {
      $result = fputcsv($fd, array($date, $remote_addr, $download_filename));
      fclose($fd);
   
      if($result > 0)
        return array(status => true);  
      else
        return array(status => false, message => 'Unable to write to '.$log_filename.'!');
    }
    else {
      return array(status => false, message => 'Unable to open log '.$log_filename.'!');
    }
  }

  date_default_timezone_set('UTC');
  # Set names of files available for download
   if (isset($_POST['release_button'])) {
     $FileNameBaseArray = array(
       "Source" => "crux-3.2.Source.tar.gz",
       "Linux32" => "crux-3.2.Linux.i686.zip",
       "Linux64" => "crux-3.2.Linux.x86_64.zip",
       "OSX" => "crux-3.2.Darwin.x86_64.zip",
       "Windows64" => "crux-3.2.Windows.AMD64.zip",
       "Windows32" => "crux-3.2.Windows.i386.zip",
     );
     $directory = "crux-3.2/";
   }
    else {
      # File names for daily build include abbreveiated unique git commit number.
      # Read commit number from file.
      $version = chop(file_get_contents("daily/latest-build.txt"));
      $FileNameBaseArray = array(
        "Source" => "crux-3.2.$version.Source.tar.gz",
        "Linux32" => "crux-3.2.$version.Linux.i686.zip",
        "Linux64" => "crux-3.2.$version.Linux.x86_64.zip",
        "OSX" => "crux-3.2.$version.Darwin.x86_64.zip",
        "Windows64" => "crux-3.2.$version.Windows.AMD64.zip",
        "Windows32" => "crux-3.2.$version.Windows.i386.zip",
      );
      $directory = "daily/";
    }
    # Redriect to file to be downloaded.
	  $downloadType = $_POST['downloadtype'];
	  $filename = $FileNameBaseArray[$downloadType];
    $filepath = $directory . $filename;
    log_download($filepath);
    header("Location: http://noble.gs.washington.edu/crux-downloads/" . $filepath);
    exit;
 ?> 
