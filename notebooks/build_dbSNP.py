import pandas as pd
import sqlite3
from tqdm import tqdm
import pysam
import os


if __name__ == "__main__":
    # base_parts = ['C:', 'Users', 'gqu', 'OneDrive - UTHealth Houston', 'projects', 'Genevic']
    #
    # # Convert path for the current operating system
    # base_path = os.path.join(*base_parts)
    # data_path = os.path.join(base_path, 'data')



    # ----------------------------
    # Configuration and Setup
    # ----------------------------

    # Paths to the VCF file and the SQLite database

    # ----------------------------
    # Configuration and Setup
    # ----------------------------

    # Paths to the VCF file and the SQLite database
    vcf_path = '../data/All_20180418.filtered.vcf.gz'  # Use the compressed VCF file
    db_path = '../data/variants.db'

    # Batch size for database inserts
    batch_size = 100000  # Adjust based on your system's capacity

    # Ensure the VCF file and its index exist
    if not os.path.exists(vcf_path):
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")

    if not os.path.exists(vcf_path + '.tbi'):
        raise FileNotFoundError(f"Index file not found: {vcf_path}.tbi")

    # ----------------------------
    # Database Setup
    # ----------------------------

    # Connect to the SQLite database (it will be created if it doesn't exist)
    conn = sqlite3.connect(db_path)
    c = conn.cursor()

    # Create table for variants with primary key on chrom and pos to prevent duplicates
    c.execute('''
        CREATE TABLE IF NOT EXISTS variants (
            chrom TEXT,
            pos INTEGER,
            id TEXT,
            ref TEXT,
            PRIMARY KEY (chrom, pos)
        )
    ''')

    # Create table to keep track of processed chromosomes and positions
    c.execute('''
        CREATE TABLE IF NOT EXISTS processed_chromosomes (
            chrom TEXT PRIMARY KEY,
            last_pos INTEGER
        )
    ''')

    # ----------------------------
    # VCF File Setup
    # ----------------------------

    # Open the compressed and indexed VCF file using pysam
    vcf = pysam.VariantFile(vcf_path)

    # Retrieve the list of chromosomes from the VCF header
    chromosomes = list(vcf.header.contigs)

    # ----------------------------
    # Processing Logic
    # ----------------------------

    # Iterate over each chromosome
    for chrom in chromosomes:
        # Check if the chromosome has been processed before
        c.execute('SELECT last_pos FROM processed_chromosomes WHERE chrom = ?', (chrom,))
        result = c.fetchone()
        if result:
            last_pos = result[0]
            print(f"Resuming chromosome {chrom} from position {last_pos}")
        else:
            last_pos = 0
            print(f"Processing chromosome {chrom} from the beginning")

        try:
            # Start a transaction for the current chromosome
            conn.execute('BEGIN')
            batch_data = []

            # Fetch records starting from the last processed position
            try:
                records = vcf.fetch(chrom, start=last_pos)
            except Exception as e:
                print(f"Error fetching records for chromosome {chrom}: {e}")
                continue  # Skip to the next chromosome

            for rec in tqdm(records, desc=f"Inserting records for {chrom}"):
                try:
                    # Append the record data to the batch
                    batch_data.append((rec.chrom, rec.pos, rec.id, rec.ref))
                    last_pos = rec.pos  # Update the last processed position
                except Exception as e:
                    print(f"Error processing record at {chrom}:{last_pos}: {e}")
                    continue  # Skip this record and continue

                # Insert batch into the database when batch_size is reached
                if len(batch_data) >= batch_size:
                    c.executemany('''
                        INSERT OR IGNORE INTO variants (chrom, pos, id, ref)
                        VALUES (?, ?, ?, ?)
                    ''', batch_data)
                    # Update last_pos in the tracking table
                    c.execute('''
                        INSERT OR REPLACE INTO processed_chromosomes (chrom, last_pos)
                        VALUES (?, ?)
                    ''', (chrom, last_pos))
                    conn.commit()
                    batch_data = []  # Reset the batch

            # Insert any remaining records after finishing the chromosome
            if batch_data:
                c.executemany('''
                    INSERT OR IGNORE INTO variants (chrom, pos, id, ref)
                    VALUES (?, ?, ?, ?)
                ''', batch_data)
                c.execute('''
                    INSERT OR REPLACE INTO processed_chromosomes (chrom, last_pos)
                    VALUES (?, ?)
                ''', (chrom, last_pos))
                conn.commit()

            # Mark chromosome as fully processed
            c.execute('''
                INSERT OR REPLACE INTO processed_chromosomes (chrom, last_pos)
                VALUES (?, ?)
            ''', (chrom, last_pos))
            conn.commit()

        except Exception as e:
            conn.rollback()  # Roll back any changes made during this chromosome
            print(f"Error processing chromosome {chrom} at position {last_pos}: {e}")
            # Update last_pos even if an error occurs
            c.execute('''
                INSERT OR REPLACE INTO processed_chromosomes (chrom, last_pos)
                VALUES (?, ?)
            ''', (chrom, last_pos))
            conn.commit()
            continue  # Proceed to the next chromosome

    # ----------------------------
    # Cleanup
    # ----------------------------

    # Close the VCF file
    vcf.close()

    # Close the database connection
    conn.close()

    print("Processing completed.")

#
#
# def initialize_database(db_path):
#     """ Initialize the SQLite database connection. """
#     return sqlite3.connect(db_path)
#
#
#
# def load_data_to_db(csv_file, conn, chunksize=500000):
#     """ Load data from CSV to SQLite database using chunks with a progress bar. """
#     # Specify delimiter as '\t' for tab-delimited files
#     use_columns = ["#chrom", "chromStart", "chromEnd", "name"]
#     total_rows = sum(1 for row in open(csv_file, 'r', encoding='utf-8')) - 1
#
#     # Create a progress bar
#     progress_bar = tqdm(total=total_rows, desc='Loading Data', unit='rows')
#
#     # Read and write data in chunks
#     for chunk in pd.read_csv(csv_file, delimiter='\t', chunksize=chunksize, usecols=use_columns):
#         chunk.to_sql('genome_data', conn, if_exists='append', index=False)
#         progress_bar.update(len(chunk))
#
#     # Close the progress bar
#     progress_bar.close()
#
# def create_indexes(conn):
#     """ Create necessary indexes on the table for faster queries. """
#     cursor = conn.cursor()
#     cursor.execute('''
#     CREATE INDEX IF NOT EXISTS idx_lookup ON genome_data ("#chrom", chromStart, chromEnd);
#     ''')
#     conn.commit()
#
#
# def batch_query(conn, query_data):
#
#     """Perform batch queries to retrieve names based on chromosomal data."""
#     cursor = conn.cursor()
#     results = []
#     query = '''
#     SELECT name FROM genome_data
#     WHERE "#chrom" = ? AND chromStart = ? AND chromEnd = ?
#     '''
#     try:
#         for params in query_data:
#             print(f"Executing query with parameters: {params}")  # Debugging output
#             cursor.execute(query, params)
#             result = cursor.fetchone()  # Fetch result for this parameter set
#             print(f"Result for {params}: {result}")  # Debugging output
#             if result:
#                 results.append(result[0])
#             else:
#                 results.append(0)
#
#     except Exception as e:
#         print(f"Error during batch query execution: {e}")
#     return results
#
#
# def close_connection(conn):
#     """ Close the database connection. """
#     conn.close()
#
#
# def setup_database(csv_file, db_path=r"C:\Users\gqu\OneDrive - UTHealth Houston\projects\Genevic\data\genomic_data.db"):
#     """ Function to build the database with data from the CSV file. """
#     conn = initialize_database(db_path)
#     load_data_to_db(csv_file, conn)
#     create_indexes(conn)
#     close_connection(conn)
#
#
# def use_database(query_data, db_path=r"C:\Users\gqu\OneDrive - UTHealth Houston\projects\Genevic\data\genomic_data.db"):
#     """ Function to use the existing database for queries. """
#     conn = initialize_database(db_path)
#
#     names = batch_query(conn, query_data)
#     print(f'Retrieved Names: {names}')
#
#     close_connection(conn)
#
# def check_csv_headers(csv_file):
#     with open(csv_file, 'r') as file:
#         first_line = file.readline()
#         print("Headers in the file:", first_line.strip())
#
#
# def count_csv_rows(csv_file):
#     """ Count the actual number of rows in the CSV file, excluding the header. """
#     with open(csv_file, 'r', encoding='utf-8') as file:
#         row_count = sum(1 for row in file) - 1  # Subtract 1 for the header
#     return row_count
#
# # Main execution function
# def main(build_db=False,db_path=r"C:\Users\gqu\OneDrive - UTHealth Houston\projects\Genevic\data\genomic_data.db"):
#     csv_file = r'C:\Users\gqu\OneDrive - UTHealth Houston\projects\Genevic\data\dbSNP_155_hg38.csv'
#     heck_csv_headers(csv_file)
#     row_count = count_csv_rows(csv_file)
#     print('number of rows = ', row_count )
#     if build_db:
#         print("Setting up the database with new data...")
#         setup_database(csv_file)
#
#
#
# if __name__ == '__main__':
#     # Set build_db to True if you need to rebuild the database with new data
#     # main(build_db=True)
#     use_database(query_data=[('chr4', 10266316, 10266317
#                               )])
