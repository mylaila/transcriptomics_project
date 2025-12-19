# Configuration VS Code pour transcriptomics_project

## Configuration automatique de l'environnement

Ce projet est configuré pour utiliser automatiquement l'environnement conda **transcriptomics**.

### Ce qui est configuré automatiquement :

1. **Interpréteur Python** : L'environnement `transcriptomics` est défini comme interpréteur par défaut
2. **Kernel Jupyter** : Le kernel "Python (transcriptomics)" est automatiquement sélectionné pour les notebooks
3. **Terminal intégré** : L'environnement `transcriptomics` est activé automatiquement à l'ouverture d'un nouveau terminal

### Vérification

Pour vérifier que tout fonctionne correctement :

1. Ouvrez un notebook → Le kernel devrait être "Python (transcriptomics)"
2. Ouvrez un nouveau terminal → Vous devriez voir `(transcriptomics)` dans le prompt
3. Dans un terminal, tapez `python --version` pour vérifier la version Python de l'environnement

### Réinstallation du kernel (si nécessaire)

Si le kernel n'apparaît pas :

```powershell
conda activate transcriptomics
python -m ipykernel install --user --name transcriptomics --display-name "Python (transcriptomics)"
```

Puis rechargez VS Code (Ctrl+R ou F1 → "Reload Window").

### Fichiers de configuration

- `settings.json` : Configuration VS Code pour Python et Jupyter
- `activate_env.ps1` : Script PowerShell d'activation (si nécessaire)
- `notebook_template.ipynb` : Template de notebook avec le bon kernel

