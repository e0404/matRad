	/* Stock GA code */
	(function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
	(i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
	m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
	})(window,document,'script','https://www.google-analytics.com/analytics.js','ga');

	ga('create', 'UA-64247123-1', 'auto');
	
	/* send custom dimensions to GA */
	function sendCustomDims(actionString) {
	
		/* Set GA clientId and print to screen */
		ga(function(tracker) {
			tracker.set('dimension1',tracker.get('clientId'))
		});
		ga(function(tracker) {
			console.log("clientId = " + tracker.get('clientId'));
		});
				
		/* Assign a timestamp to the utc_millisecs custom dimension */
		timeTag = new Date().getTime();
		ga('set', 'dimension2', timeTag);
		console.log("time = " + timeTag + " ms");
	
		/* identify action */
		ga('set', 'dimension3',actionString)
		console.log('action = ' + actionString);
		
		/* Send a pageview */
		ga('send', 'pageview');
	
		return true;
	}
	
	sendCustomDims('landing');		
	
	console.log("ga init done");
	
	// implementation according to https://stackoverflow.com/questions/688196/how-to-use-a-link-to-call-javascript
	// Wait for the page to load first
	window.onload = function() {

	  //Get a reference to the link on the page
	  // with an id of "mylink"
	  var pdfoverview     = document.getElementById("pdfoverview");
	  var downloadZipball = document.getElementById("downloadZipball");
	  var downloadTarball = document.getElementById("downloadTarball");
	  var downloadStandalone = document.getElementById("downloadStandalone");
	  var brownbagYoutube = document.getElementById("brownbagYoutube");
	  var gotoGithub      = document.getElementById("gotoGithub");
	  var gotoMap         = document.getElementById("gotoMap");
	  var gotoWiki        = document.getElementById("gotoWiki");
	  var gotoMatlab      = document.getElementById("gotoMatlab");
	  var gotoDKFZ        = document.getElementById("gotoDKFZ");
	  var gotoE0404       = document.getElementById("gotoE0404");
	  var gotoAuthorlist  = document.getElementById("gotoAuthorlist");
	  var gotoLicense     = document.getElementById("gotoLicense");
	  var gotoIpopt       = document.getElementById("gotoIpopt");
	  var writeEmail      = document.getElementById("writeEmail");

	  pdfoverview.onclick     = function() {
												sendCustomDims('pdfoverview');
												return true;
											}
	  downloadZipball.onclick = function() {
												sendCustomDims('downloadZipball');
												return true;
											}
	  downloadTarball.onclick = function() {
												sendCustomDims('downloadTarball');
												return true;
											}
	  downloadStandalone.onclick = function() {
												sendCustomDims('downloadStandalone');
												return true;
											}
	  brownbagYoutube.onclick = function() {
												/* xxx this does not trigger */
												return true;
											}
	  gotoGithub.onclick      = function() {
												sendCustomDims('gotoGithub');
												return true;
											}
	  gotoMap.onclick         = function() {
												sendCustomDims('gotoMap');
												return true;
											}
	  gotoWiki.onclick        = function() {
												sendCustomDims('gotoWiki');
												return true;
											}
	  gotoMatlab.onclick      = function() {
												sendCustomDims('gotoMatlab');
												return true;
											}
	  gotoDKFZ.onclick        = function() {
												sendCustomDims('gotoDKFZ');
												return true;
											}
	  gotoE0404.onclick       = function() {
												sendCustomDims('gotoE0404');
												return true;
											}
	  gotoAuthorlist.onclick  = function() {
												sendCustomDims('gotoAuthorlist');
												return true;
											}
	  gotoLicense.onclick     = function() {
												sendCustomDims('gotoLicense');
												return true;
											}
	  gotoIpopt.onclick       = function() {
												sendCustomDims('gotoIpopt');
												return true;
											}
	  writeEmail.onclick      = function() {
												sendCustomDims('writeEmail');
												return true;
											}
	}