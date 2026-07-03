# GitHub Desktop steps

## Current state

If GitHub Desktop shows `No local changes`, your files have already been committed locally.
That is normal.

## Recommended workflow

1. Run the project locally.

```bash
python main.py
python -m pytest
```

2. Check generated figures in:

```text
reports/figures/
```

3. Open GitHub Desktop.

4. If there are changes, write a summary such as:

```text
Add R-style transformed-axis reliability plots
```

5. Click:

```text
Commit to main
```

6. When you are ready to put the project online, click:

```text
Publish repository
```

7. Keep the repository public if you plan to put it on your resume.

8. After publishing, future updates use:

```text
Commit to main -> Push origin
```
