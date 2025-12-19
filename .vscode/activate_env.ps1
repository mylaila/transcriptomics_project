# Script d'activation automatique de l'environnement transcriptomics
# Ce script est exécuté automatiquement lors de l'ouverture d'un terminal

Write-Host "🔄 Activation de l'environnement conda 'transcriptomics'..." -ForegroundColor Cyan
conda activate transcriptomics

if ($LASTEXITCODE -eq 0) {
    Write-Host "✅ Environnement 'transcriptomics' activé avec succès!" -ForegroundColor Green
    Write-Host "📍 Répertoire de travail: $PWD" -ForegroundColor Yellow
} else {
    Write-Host "❌ Erreur lors de l'activation de l'environnement" -ForegroundColor Red
}
