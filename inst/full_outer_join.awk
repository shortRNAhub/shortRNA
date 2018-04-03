#!/usr/bin/awk -f
# NAME
#     full-outer-join.awk - perform a full outer join on tabular files
#
# SYNOPSIS 
#     full-outer-join.awk [HEADER=1] [SORT_OPTIONS=sort_options] FILE...
#
# DESCRIPTION
#     Perform a full outer join on the first field of two or more delimited
#     (default: tab) text files.
#
# OPERANDS
#     FILE
#         A delimited (default: tab) text file with 1 header line. The first
#         field of each record will be used as the the join field, must be
#         unique within each file, and needn't be sorted.
#
#         The field separator (delimiter) can be changed by adding an
#         FS="DELIMITER" operand before FILE. The new field separator will be
#         in effect for the remaining files, or until changed.
#
#     HEADER
#         Set to "1" (or any non-zero, non-blank value) if the first row is a
#         header.
#
#     SORT_OPTIONS
#         Options to pass to the "sort" command. By default, "-k 1,1" is used
#         to sorted output lexically on the first (join) field. Specify, e.g.,
#         "-k 1n,1n" to sort output numerically on the join field.
#
# STDIN
#     Standard input will be used only if no file operands are specified, or
#     if a file operand is '-'.
#
# STDOUT
#     A full outer join on the input files, with any headers printed in the
#     order in which they appear in the input files. The join field will be
#     output only once. Null values are printed as "NULL". By default, fields
#     are separated by a tab.
#
# STDERR
#     An error message will be printed to standard error, and the program will
#     terminate with a non-zero exit status, if one of the following occurs:
#
#     1. The number of input files is 1.
#     2. A duplicate join field value in the same input file is detected.
#     3. The number of fields in a record is inconsistent.
#
# EXAMPLES
#     1. Join 3 tab-delimited text files containing a header:
#
#        $ cat file1.txt
#        join_field	column1	column2
#        XYZ	c11	c12
#        ABC	c21	c22
#        GOP	c31	c32
#        $ cat file2.txt
#        join_field	column3	column4
#        GOP	c33	c34
#        XYZ	c13	c14
#        BOG	c43	c44
#        $ cat file3.txt
#        join_field	column5	column6
#        XYZ	c15	c16
#        GOP	c35	c36
#        BOG	c45	c46
#        $ full-outer-join.awk HEADER=1 file1.txt file2.txt file3.txt
#        join_field	column1	column2	column3	column4	column5	column6
#        ABC	c21	c22	NULL	NULL	NULL	NULL
#        BOG	NULL	NULL	c43	c44	c45	c46
#        GOP	c31	c32	c33	c34	c35	c36
#        XYZ	c11	c12	c13	c14	c15	c16
#
#    2. Join 3 comma-delimited text files that have "old Mac" newlines (i.e.,
#       carriage return instead of line feed) as is exported by Excel for Mac:
#       $ full-outer-join.awk FS="," RS="\r" HEADER=1 file[123].txt
#       ...
#
#    3. Join tab-delimited files containing SNPs on the same chromosome from
#       three different samples, where the file format is
#       POSITION REFERENCE_BASE SNP_BASE
#       Output is sorted numerically by position. A shell supporting process
#       substitution (such as bash, ksh, or zsh) is assumed.
#
#       $ cat sample1.snp
#       100	A	C
#       200	T	G
#       300	G	T
#       $ cat sample2.snp
#       100	A	C
#       250	T	G
#       300	G	T
#       $ cat sample3.snp
#       100	A	C
#       275	T	G
#       300	G	A
#       $ full-outer-join.awk <(cut -f 1,2 sample[123].snp | sort -u) \
#                             <(cut -f 1,3 sample1.snp) \
#                             <(cut -f 1,3 sample2.snp) \
#                             <(cut -f 1,3 sample3.snp) | sort -k 1n,1n
#       100	A	C	C	C
#       200	T	G	NULL	NULL
#       250	T	NULL	G	NULL
#       275	T	NULL	NULL	G
#       300	G	T	T	A
#      
# CHANGE HISTORY
#     2014-03-04    Fix header-display bug where only the header from the
#                   first file was displayed.
#     2013-01-29    Bug fixes for changes in 2013-01-28 version.
#     2013-01-28    Make the header row optional, disabled by default and 
#                   enabled by the "HEADER=1" operand. Add a SORT_OPTIONS
#                   operand, which is useful when "HEADER=1" and, e.g.,
#                   the join field is to be sorted numerically instead of
#                   lexically.
#     2012-11-02    Basic error checking.
#     2012-11-01    Initial version
#
# BUGS
#     Minimal error checking is done on the input.
#
# AUTHOR
#     Nathan Weeks <nathan.weeks@ars.usda.gov>

BEGIN { FS = OFS = "\t" }

# 1st column contains primary key
NR == 1 && HEADER { header[1] = $1 }

FNR == 1 {
    num_fields_in_file = NF
    if (HEADER)
        for (column = 2; column <= NF; column++)
            header[column+(degree == 0 ? 0 : degree - 1)] = $column
    degree += (degree ? NF-1 : NF) # count join field once
}

!HEADER || (HEADER && FNR > 1) {
    # error if the number of fields differs from the number of header fields
    if (NF != num_fields_in_file) {
        print FILENAME ":" FNR \
              ": number of fields (" NF ") differs from the number of " \
              "fields in first record (" num_fields_in_file"):" \
              | "cat 1>&2"
        print | "cat 1>&2"
        exit_status = 1; exit 1
    }

    # error if we've already seen this join field value before in this file
    if (join_fields[$1] > NR - FNR) {
        print FILENAME ":" FNR ": duplicate join field: " $1 | "cat 1>&2"
        exit_status = 1; exit 1
    } else
        join_fields[$1] = NR

    for (column = 2; column <= NF; column++)
        row[$1,degree-NF+column] = $column # save remaining columns
}

END {
    if (exit_status) # if abnormal termination
        exit exit_status

    if (FNR == NR) {
        print "ERROR: only one input file; must have two or more" | "cat 1>&2"
        exit 1
    }

    if (HEADER) {
        for (column = 1; column < degree; column++)
            printf("%s\t", header[column])
        printf("%s\n", header[degree])
    }

    sort_command = "sort " (SORT_OPTIONS ? SORT_OPTIONS : "-k 1,1")

    # print data rows
    for (key in join_fields) {
        printf("%s", key) | sort_command
        for (column = 2; column <= degree; column++)
            if ((key, column) in row)
                printf("\t%s", row[key,column]) | sort_command
            else
                printf("\tNULL") | sort_command
        printf("\n") | sort_command
    }
}

