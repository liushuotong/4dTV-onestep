#ifndef DATA_READER_H
#define DATA_READER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE_LENGTH 1024
#define INITIAL_CAPACITY 100

// TSV Record Structure
typedef struct {
    char *gene;
    char *paralog_gene;
} TsvRecord;

// FASTA Record Structure
typedef struct {
    char *header;
    char *sequence;
} FastaRecord;

// Function prototypes for TSV file operations
void read_tsv(const char*, TsvRecord**, int*);
void free_tsv_records(TsvRecord*, int);
TsvRecord* tsv_reader(const char*, int*);

// Function prototypes for FASTA file operations
void read_fasta(const char*, FastaRecord**, int*);
void free_fasta_records(FastaRecord*, int);
FastaRecord* fasta_reader(const char*, int*);


/**
 * Reads a TSV file and stores the records in an array of TsvRecord structs.
 * The TSV file is expected to have two columns:
 * the first for the gene and the second for the paralog_gene.
 */
void read_tsv(const char *filename, TsvRecord **records, int *num_records) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    TsvRecord *temp_records = (TsvRecord *)malloc(INITIAL_CAPACITY * sizeof(TsvRecord));
    if (temp_records == NULL) {
        perror("Error allocating memory");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    char line[MAX_LINE_LENGTH];
    int record_count = 0;
    int record_capacity = INITIAL_CAPACITY;

    while (fgets(line, sizeof(line), file)) {
        line[strcspn(line, "\r\n")] = '\0'; // Remove newline characters

        char *gene = strtok(line, "\t");
        char *paralog_gene = strtok(NULL, "\t");

        if (gene == NULL || paralog_gene == NULL) {
            continue; // Skip lines that don't have both columns
        }

        TsvRecord current_record;
        current_record.gene = strdup(gene);
        current_record.paralog_gene = strdup(paralog_gene);

        if (current_record.gene == NULL || current_record.paralog_gene == NULL) {
            perror("Error duplicating string");
            fclose(file);
            free_tsv_records(temp_records, record_count);
            exit(EXIT_FAILURE);
        }

        temp_records[record_count++] = current_record;

        if (record_count >= record_capacity) { // Reallocate memory if needed
            record_capacity *= 2;
            TsvRecord *new_records = (TsvRecord*)realloc(temp_records, record_capacity * sizeof(TsvRecord));
            if (new_records == NULL) {
                perror("Error reallocating memory");
                fclose(file);
                free_tsv_records(temp_records, record_count);
                exit(EXIT_FAILURE);
            }
            temp_records = new_records;
        }
    }

    fclose(file);

    *records = temp_records;
    *num_records = record_count;
}

/**
 * Frees the memory allocated for a list of TSV records.
 */
void free_tsv_records(TsvRecord *records, int num_records) {
    for (int i = 0; i < num_records; i++) {
        free(records[i].gene);
        free(records[i].paralog_gene);
    }
    free(records);
}

/**
 * Reads a TSV file and returns the records without freeing them.
 */
TsvRecord* tsv_reader(const char* location, int *num_records) {
    TsvRecord *records;
    read_tsv(location, &records, num_records);
    return records;
}


/**
 * Reads a FASTA file and stores the records in an array of FastaRecord structs.
 */
void read_fasta(const char *filename, FastaRecord **records, int *num_records) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    FastaRecord *temp_records = (FastaRecord *)malloc(INITIAL_CAPACITY * sizeof(FastaRecord));
    if (temp_records == NULL) {
        perror("Error allocating memory");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    char line[MAX_LINE_LENGTH];
    int record_count = 0;
    int record_capacity = INITIAL_CAPACITY;

    FastaRecord current_record = {NULL, NULL};

    while (fgets(line, sizeof(line), file)) {
        line[strcspn(line, "\r\n")] = '\0'; // Remove newline characters

        if (line[0] == '>') { // New record starts
            if (current_record.header != NULL) { // Save previous record
                temp_records[record_count++] = current_record;

                if (record_count >= record_capacity) { // Reallocate memory if needed
                    record_capacity *= 2;
                    FastaRecord *new_records = (FastaRecord*)realloc(temp_records, record_capacity * sizeof(FastaRecord));
                    if (new_records == NULL) {
                        perror("Error reallocating memory");
                        fclose(file);
                        free_fasta_records(temp_records, record_count);
                        exit(EXIT_FAILURE);
                    }
                    temp_records = new_records;
                }
            }

            // Start new record
            current_record.header = strdup(line + 1);
            if (current_record.header == NULL) {
                perror("Error duplicating string");
                fclose(file);
                free_fasta_records(temp_records, record_count);
                exit(EXIT_FAILURE);
            }
            current_record.sequence = NULL;
        } else { // Append sequence
            if (current_record.sequence == NULL) {
                current_record.sequence = strdup(line);
            } else {
                size_t new_length = strlen(current_record.sequence) + strlen(line) + 1;
                char *new_sequence = (char *)realloc(current_record.sequence, new_length);
                if (new_sequence == NULL) {
                    perror("Error reallocating memory");
                    fclose(file);
                    free_fasta_records(temp_records, record_count);
                    free(current_record.header);
                    free(current_record.sequence);
                    exit(EXIT_FAILURE);
                }
                current_record.sequence = new_sequence;
                strcat(current_record.sequence, line);
            }

            if (current_record.sequence == NULL) {
                perror("Error duplicating string");
                fclose(file);
                free_fasta_records(temp_records, record_count);
                free(current_record.header);
                exit(EXIT_FAILURE);
            }
        }
    }

    if (current_record.header != NULL) { // Save the last record
        temp_records[record_count++] = current_record;
    }

    fclose(file);

    *records = temp_records;
    *num_records = record_count;
}

/**
 * Frees the memory allocated for a list of FASTA records.
 */
void free_fasta_records(FastaRecord *records, int num_records) {
    for (int i = 0; i < num_records; i++) {
        free(records[i].header);
        free(records[i].sequence);
    }
    free(records);
}

/**
 * Reads a FASTA file and returns the records without freeing them.
 */
FastaRecord* fasta_reader(const char* location, int *num_records) {
    FastaRecord *records;
    read_fasta(location, &records, num_records);
    return records;
}

#endif // DATA_READER_H
