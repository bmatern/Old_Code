SET CodePath=C:\MinIONScripts\FixSwisslabFiles
SET SpecFile=SwisslabInstallerOptions_Windows.spec
SET CondaEnvironment=AlleleSubEnvironment

:: Run Pyinstaller to create executables
cd %CodePath%
activate %CondaEnvironment% && pyinstaller %SpecFile% && deactivate
