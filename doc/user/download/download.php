<?php

# <OPTION value="Linux32">&nbsp;&nbsp; Linux 32 bit binaries (NO vendor reader support)</OPTION>
# <OPTION value="Linux64">&nbsp;&nbsp; Linux 64 bit binaries (NO vendor reader support)</OPTION>
# <OPTION value="OSX">&nbsp;&nbsp; Macintosh binaries (NO vendor reader support)</OPTION>
# <OPTION value="Windows32">&nbsp;&nbsp; Windows (includes vendor reader support)</OPTION>
# <OPTION value="Source">&nbsp;&nbsp; Source, CMake build (includes vendor reader support)</OPTION>

$debug = 0;

$FileNameBaseArray = array(
	"Source" => "crux-2.0/crux-2.0.tar.bz2", #pwiz-src-3_0_3725.tar.bz2 
	"Linux32" => "crux-2.0/crux-2.0.Linux.i686.zip", #pwiz-src-without-v-3_0_3725.tar.bz2  
	"Linux64" => "crux-2.0/crux-2.0.Linux.x86_64.zip", # libpwiz_msvc_3_0_3923.zip
	"OSX" => "crux-2.0/crux-2.0-.Darwin.x86_64.zip", # libpwiz_src_3_0_3923.tgz
);


#echo "hello<BR>\n";
$count_my_page = ("hitcounter_plus.txt");
#$hits = file($count_my_page);
#$counts = $hits[0];
#$counts+=1;
$fp = fopen($count_my_page,"a") or die ("Error opening file in write mode!<BR>");

$today = date("F j, Y, g:i a");
#echo $today," hello<BR>\n";

######$downloadType = $_GET['downloadtype'];
$downloadType = $_POST['downloadtype'];

$putString =  $today."\t".$_SERVER['REMOTE_HOST']."\t".$_SERVER['REMOTE_ADDR']."\t".$downloadType."\n";

#echo "hello<BR>\n";

#echo $putString,"<BR>\n";

fputs($fp , $putString);
fclose($fp);

#  file_get_contents("http://proteowizard.sourceforge.net/downloadstub.php");

if ($downloadType) { // if page is not submitted to itself echo the form
	#echo "trying to download\n";

	$winInstaller = 0;	

	$downloadTypeString = $downloadType;

	if(($downloadType == 'bt36i') || ($downloadType == 'bt83i')){
		$winInstaller = 'i';
		$downloadTypeString = rtrim($downloadType,"i");
	}
	
	if($downloadType == 'bt81n'){
		$downloadTypeString = 'bt81';
	}

	if($debug == 1){
		echo "DOWNLOAD TYPE STRING: ",$downloadTypeString,"<BR>\n";
	}

    $remoteURL = "http://teamcity.labkey.org:8080/app/rest/buildTypes/id:" . $downloadTypeString . "/builds?status=SUCCESS&count=1&guest=1";

//    $remoteURL = "http://teamcity.labkey.org:8080/app/rest/buildTypes/id:bt36/builds?status=SUCCESS&count=1&guest=1";

	if($debug == 1){
		echo "<BR>REMOTE URL: ".$remoteURL."<BR>\n";
	}

	$teamCityInfoString = file_get_contents("$remoteURL");
		 
	preg_match("/build id=\"(\d+)\"/", $teamCityInfoString, $matches);
	$buildId = $matches[1];
	$versionURL = "http://teamcity.labkey.org:8080/repository/download/" . $downloadTypeString . "/" . $buildId . ":id/VERSION?guest=1";
	
	if($debug == 1){
		echo "VERSION URL: ",$versionURL,"<BR>";
	}
	
	$versionString = file_get_contents($versionURL);
	if(!$winInstaller){
		$versionString = preg_replace('/\./', '_', $versionString);
		$downloadURL = "http://teamcity.labkey.org:8080/repository/download/" . $downloadTypeString . "/". $buildId . ":id/" . $FileNameBaseArray[$downloadType];
		
	}
	else{
		$downloadURL = "http://teamcity.labkey.org:8080/repository/download/" . $downloadTypeString . "/". $buildId . ":id/" . $FileNameBaseArray[$downloadType];
	
	}
	
	$downloadURL = preg_replace("/XXX/", $versionString, $downloadURL);
	$downloadFile = basename($downloadURL);
	$downloadURL = $downloadURL."?guest=1";

	if($debug == 1){
		echo "DOWNLOAD URL: ",$downloadURL,"<BR>\n";
	}
	
#	$artifactURL = "http://teamcity.labkey.org:8080/viewLog.html?buildId=". $buildId ."&tab=artifacts&buildTypeId=". $downloadType . "&guest=1";
#	$artifactString = file_get_contents($artifactURL);
#	$matchString = '/<a href=\'repository\/download\/' . $downloadType . '\/' . $buildId . ':id\/(.*?)\'>/';
#	preg_match($matchString,$artifactString, $matches);
#	$downloadFile = $matches[1];

#	echo "$artifactURL<BR>\n";

#	preg_match_all($matchString,$artifactString, $matches);
#	if($winInstaller == 1){
#	$downloadFile = $matches[1][1];
#	}
#	else{
#		$downloadFile = $matches[1][0];
#	}


#	if($noVendor == 1){
#		$downloadFile = str_replace('pwiz-src-','pwiz-src-without-v-',$downloadFile);
#	}
	
	
	
#	$downloadURL = "http://teamcity.labkey.org:8080/guestAuth/repository/download/" . $downloadType . "/" .  	$buildId . ":id/" . $downloadFile ;
#	echo "DOWNLOAD URL:".$downloadURL."<BR>";
if($debug == 1){
	exit();
}

    $mm_type="application/x-bzip2"; 

	header("Pragma: public");
	header("Expires: 0");
	header("Cache-Control: must-revalidate, post-check=0, pre-check=0");
	header("Cache-Control: public");
	header("Content-Description: File Transfer");
	header("Content-Type: " . $mm_type);
#	header("Content-Length: " .(string)($fs[1]) );  //for whatever reason I can't the size...
	header('Content-Disposition: attachment; filename="'.$downloadFile.'"');
	header("Content-Transfer-Encoding: binary\n");
	readfile($downloadURL); // outputs the content of the file
	exit;
}
else {
	print 'Please go <a href="http://proteowizard.sourceforge.net/downloads.shtml">HERE</a> to download ProteoWizard' ;
}

?>
