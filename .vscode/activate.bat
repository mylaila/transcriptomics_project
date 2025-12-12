@echo off
REM Script pour activer l'environnement conda transcriptomics
echo Activation de l'environnement conda 'transcriptomics'...
call conda activate transcriptomics
if %errorlevel% equ 0 (
    echo [OK] Environnement 'transcriptomics' active avec succes!
    echo Repertoire: %CD%
) else (
    echo [ERREUR] Echec de l'activation de l'environnement
)
