.PHONY: run-notebooks

run-notebooks:
	jupyter lab --allow-root --no-browser --port 8888 --NotebookApp.custom_display_url=http://${CODESPACE_NAME}-8888.app.github.dev