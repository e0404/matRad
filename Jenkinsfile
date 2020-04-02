node {
    def osName
    stage('Set Environment') {
        echo 'Setting the environment..'
		checkout scm
        if (isUnix()) {
            def uname = sh script: 'uname', returnStdout: true
            if (uname.startsWith("Darwin")) {
                osName= "Macos"
				sh label: '', script: 'chmod +x ./submodules/MCsquare/MCsquare_mac'
                sh label: '', script: 'chmod +x runtests.sh'
            }
            else {
                osName= "Linux"
                sh label: '', script: 'chmod +x ./submodules/MCsquare/MCsquare_linux'
                sh label: '', script: 'chmod +x runtests.sh'
            }
        }
        else {
            osName= "Windows"
        }

    }
    stage('Test with MATLAB') 
	{
        echo 'Testing on MATLAB..'
        echo 'Make sure you can run matlab script from the command line'
        if (osName == "Windows") {
            bat label: '', script: 'matlab -noFigureWindow -nosplash -minimize -wait -r "cd unitTest; matRad_runTests"'
        }
        else {
            sh label: '', script: 'matlab -nodisplay -r "cd unitTest; matRad_runTests"'
        } 
    }
	
	stage('Test with Octave') 
	{
		echo 'Testing on Octave..'
        if (osName == "Windows") {
            bat label: '', script: './runtests.sh octave-cli'
        }
        else {
            sh label: '', script: './runtests.sh octave-cli'
        } 
	}
	
    stage('Build') {
        echo 'Building....'
        if (osName == "Windows") {
            bat label: '', script: 'matlab -nodisplay -r "matRad_compileStandalone"'
        }
        else {
            sh label: '', script: 'matlab -nodisplay -r "matRad_compileStandalone"'
        } 
    }
    stage('Archive') {
        echo 'Archiving....'
        archiveArtifacts artifacts: 'standalone/for_redistribution/matRad_installer*', onlyIfSuccessful: true
    }
}