.PHONY: test-chm13 test install clean

# Test dataset paths
CHM13_FASTA = ~/Dropbox/workspace/new/biology/satellome/test_dataset/chm13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna
CHM13_OUTPUT = ~/Dropbox/workspace/new/biology/satellome/test_dataset/chm13/output
THREADS = 8

# Run satellome on CHM13 test dataset
test-chm13:
	@echo "Running satellome on CHM13 dataset..."
	@mkdir -p $(CHM13_OUTPUT)
	satellome -i $(CHM13_FASTA) -o $(CHM13_OUTPUT) -p chm13 -t $(THREADS)

# Run with force flag
test-chm13-force:
	@echo "Running satellome on CHM13 dataset (force rerun)..."
	@mkdir -p $(CHM13_OUTPUT)
	satellome -i $(CHM13_FASTA) -o $(CHM13_OUTPUT) -p chm13 -t $(THREADS) --force

# Run with TRF enabled
test-chm13-trf:
	@echo "Running satellome on CHM13 dataset with TRF..."
	@mkdir -p $(CHM13_OUTPUT)
	satellome -i $(CHM13_FASTA) -o $(CHM13_OUTPUT) -p chm13 -t $(THREADS) --run-trf

# Install package in development mode
install:
	pip install -e .

# Run unit tests
test:
	python -m pytest tests/ -v

# Clean build artifacts
clean:
	rm -rf build/ dist/ *.egg-info src/*.egg-info
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete 2>/dev/null || true
