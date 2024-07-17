import gzip

def read_gzip_file(path):
    with gzip.open(path, 'rt') as f:  # 'rt' mode for reading text
        contents = f.read()
        print(contents)

# Usage
read_gzip_file('/Users/katherinetian/Downloads/GWAS/covariance.txt.gz')
