#ifndef STRING_UTILS_H
#define STRING_UTILS_H

#include <ctype.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* **************************************
 *  		Functions		*
 * *************************************/

/**
*  @brief Compares two given strings
*  @param str1 string 1 to compare
*  @param str2 string 2 to compare
*  @return 0 if equal string, more than 0 if str1 is greater than str2 and 
*          less than 0 if str1 is less than str2
*  
*  Compares lexicographically two given strings
*/

int equals(const char *str1, const char *str2);

/**
*  @brief Compares two given strings not considering upper and lowercase differences
*  @param str1 string 1 to compare
*  @param str2 string 2 to compare
*  @return 0 if equal string, more than 0 if str1 is greater than str2 and 
*          less than 0 if str1 is less than str2
*  
*  Compares lexicographically two given strings not considering upper and lowercase differences 
*/
int equals_ignore_case(const char *str1, const char *str2);

/**
*  @brief Determines if a given string is numeric
*  @param str string to process
*  @return 1 if numeric, 0 otherwise
*  
*  Determines if a given string is numeric
*/
int is_numeric(const char *str);

/**
*  @brief Determines if a string fragment is present in a given string
*  @param str string to process
*  @param search string to search
*  @return 0 if not found, 1 if search found
*  
*  Determines if a string fragment is present in a given string
*/
int contains(const char *str, const char *search);

/**
*  @brief Determines if a string starts with a given string fragment
*  @param str string to process
*  @param search string to search at the start of the string
*  @return 0 string does not start with search string, 1 string starts with search string
*  
*  Determines if a string starts with a given string fragment
*/
int starts_with(const char *str, const char *search);

/**
*  @brief Determines if a string starts with a given string fragment
*  @param str string to process
*  @param search string to search at the start of the string
*  @param length maximum number of characters to check
*  @return 0 string does not start with search string, 1 string starts with search string
*  
*  Determines if a string starts with a given string fragment
*/
int starts_with_n(const char *str, const char *search, int length);

/**
*  @brief Determines if a string ends with a given string fragment
*  @param str string to process
*  @param search string to search at the end of the string
*  @return 0 string does not end with search string, 1 string ends with search string
*  
*  Determines if a string ends with a given string fragment
*/
int ends_with(const char *str, const char *search);

/**
*  @brief Converts a given string to lowercase 
*  @param str string to convert
*  @return lowercase string
*  
*  Converts a given string to lowercase 
*/
char* to_lower_case(char *str);

/**
*  @brief Converts a given string to uppercase 
*  @param str string to convert
*  @return uppercase string
*  
*  Converts a given string to uppercase 
*/
char* to_upper_case(char *str);

/**
*  @brief Sets the end of the string in one character less than actual string length
*  @param str string to process
*  @return char at the last string position
*  
*  Sets the end of the string in one character less than actual string length
*/
char chop(char *str);

/**
*  @brief Sets the end of the string in the given position
*  @param str string to process
*  @param position end position of the processed string
*  @return char at given position
*  
*  Sets the end of the string in the given position
*/
char chop_at(char *str, int position);

/**
*  @brief Substitutes '\n' an '\r' characters by '\0' in the last position of a string
*  @param str string to process
*  @return '\0' if no substitution, '\n' or '\r' replaced character if substituted
*  
*  Substitutes '\n' and '\r' characters by '\0' in the last position of a string
*/
char chomp(char *str);

/**
*  @brief Substitutes '\n' and '\r' characters by '\0' in the given position
*  @param str string to process
*  @param position substitution position
*  @return '\0' if no substitution, '\n' or '\r' replaced character if substitution is made
*  
*  Substitutes '\n' and '\r' characters by '\0' in the given position
*/
char chomp_at(char *str, int position);

/**
*  @brief Removes a given char from a string
*  @param str string to process
*  @param c char to remove from string
*  @return pointer to the processed string
*  
*  Removes a given char from a string
*/
char* remove_char(char *str, char c);

/**
*  @brief Removes the char from a string at a given position
*  @param str string to process
*  @param position char position to remove from string
*  @return pointer to the processed string
*  
*  Removes the char from a string at a given position
*/
char* remove_char_at(char *str, int position);

/**
*  @brief Removes a given string fragment from a string
*  @param str string to process
*  @param search_str string fragment to remove
*  @return Removes a given string fragment from a string
*  
*  Removes the char from a string at a given position
*/
char* remove_str(char *str, const char *search_str);

/**
*  @brief Removes a given number of characters from the start of a string
*  @param str string to process
*  @param length number of characters to remove from start of the string
*  @return pointer to the processed string
*  
*  Removes a given number of characters from the start of a string
*/
char* remove_start(char *str, int num_chars);

/**
*  @brief Removes a given number of characters from the end of a string
*  @param str string to process
*  @param length number of characters to remove from end of the string
*  @return pointer to the processed string
*  
*  Removes a given number of characters from the end of a string
*/
char* remove_end(char *str, int num_chars);

/**
*  @brief Replaces a given character by another in a string
*  @param str string to process
*  @param orig char to replace
*  @param repl new char to place on
*  @param max_line_length maximum number of characters of the string to process
*  @return pointer to the processed string
*  
*  Replaces a given character by another in a string
*/
char* str_replace(char *str, const char *orig, const char *repl, int max_line_length);

/**
*  @brief Concatenates two string arrays in one array
*  @param dest destination array
*  @param orig1_length length of the first array to concatenate
*  @param orig1 first array to concatenate
*  @param orig2_length length of the second array to concatenate
*  @param orig2 second array to concatenate
*  @return number of concatenated strings (orig1_length + orig2_length)
*  
*  Concatenates two string arrays in one array
*/
int array_concat(char **dest, int orig1_length, const char **orig1, int orig2_length, const char **orig2);

/**
*  @brief Removes spaces from left and right side of the string
*  @param str string to trim
*  @return processed string
*  
*  Removes spaces from left and right side of the strings
*  Special characters that are considered spaces are:
*      ' '    (0x20)    space (SPC)
*      '\t'   (0x09)    horizontal tab (TAB)
*      '\n'   (0x0a)    newline (LF)
*      '\v'   (0x0b)    vertical tab (VT)
*      '\f'   (0x0c)    feed (FF)
*      '\r'   (0x0d)    carriage return (CR)
*/
char* trim(char *str);

/**
*  @brief Removes spaces from left side of the string
*  @param str string to trim
*  @return processed string
*  
*  Removes spaces from left side of the strings
*  Special characters that are considered spaces are:
*      ' '    (0x20)    space (SPC)
*      '\t'   (0x09)    horizontal tab (TAB)
*      '\n'   (0x0a)    newline (LF)
*      '\v'   (0x0b)    vertical tab (VT)
*      '\f'   (0x0c)    feed (FF)
*      '\r'   (0x0d)    carriage return (CR)
*/
char* ltrim2(char *str);

/**
*  @brief Removes spaces from right side of the string
*  @param str string to trim
*  @return processed string
*  
*  Removes spaces from right side of the strings
*  Special characters that are considered spaces are:
*      ' '    (0x20)    space (SPC)
*      '\t'   (0x09)    horizontal tab (TAB)
*      '\n'   (0x0a)    newline (LF)
*      '\v'   (0x0b)    vertical tab (VT)
*      '\f'   (0x0c)    feed (FF)
*      '\r'   (0x0d)    carriage return (CR)
*/
char* rtrim2(char *str);

/**
*  @brief Cuts from left side of a string the given number or characters
*  @param string string to trim
*  @param num_chars number of characters to trim
*  @return processed string
*  
*  Cuts from left side of a string the given number or characters
*/
char* ltrim(char* string, int num_chars);

/**
*  @brief Cuts from right side of a string the given number or characters
*  @param string string to trim
*  @param num_chars number of characters to trim
*  @return processed string
*  
*  Cuts from right side of a string the given number or characters
*/
char* rtrim(char* string, int num_chars);

/**
*  @brief Removes blank spaces from left and right side of the string
*  @param str string to strip
*  @return processed string
*  
*  Removes blank spaces (' ') from left and right side of the strings 
*/
char* strip(char *str);

/**
*  @brief Removes blank spaces from left side of the string
*  @param str string to strip
*  @return processed string
*  
*  Removes blank spaces (' ') from left side of the strings
*/
char* lstrip(char *str);

/**
*  @brief Removes blank spaces from right side of the string
*  @param str string to strip
*  @return processed string
*  
*  Removes blank spaces (' ') from right side of the strings
*/
char* rstrip(char *str);

/**
*  @brief Splits a string into string fragments using a blank space separator (' ')
*  @param str string to split
*  @return pointers to the splitted strings
*  
*  Splits a string into string fragments using a blank space separator (' ')
*/
char** split(char *str, const char *delimiters, int *num_substrings);

/**
*  @brief Splits a string into string fragments using a blank space separator (' ')
*  @param str string to split
*  @param limit maximum number of splitted strings
*  @return pointers to the splitted strings
*  
*  Splits a string into string fragments using a blank space separator (' ')
*  No more than limit splitted fragments will be returned
*/
char** splitn(char *str, const char *delimiters, int limit, int *num_substrings);

/**
 */
unsigned int get_to_first_blank(char *str_p, unsigned int length, char *str_out_p);

/**
 * @brief Case-insensitive string comparison. Inspired in non-standard function:
 *        http://www.opensource.apple.com/source/gpatch/gpatch-2/patch/strcasecmp.c
 * 
 * @param s1 First string to compare
 * @param s2 Second string to compare  
 * @return Less, equals or greater than zero if s1 is lexicographically less 
 * 	   than, equal to or greater than s2.
 */
int strcasecmp(const char *s1, const char *s2);

/* ******************************************************************************
 *    		Functions to encode/decode nucleotid sequences			*
 * *****************************************************************************/

//enum bases{ DD = -1, AA = 0, CC = 1, GG = 2, TT = 3 }; //Lexicographic order
static const char alph_rep[] ={'A', 'C', 'G', 'T'};

/**
 *  @brief Inits table for nucleotide coding/decoding 
 *  @return void
 * 
 *  Inits table[128] for nucleotide coding/decoding 
 */
void initTable();

/**
 *  @brief Encodes a sequence of plain nucleotides
 *  @param dest pointer to destination char with encoded nucleotides
 *  @param src pointer to char with plain nucleotides
 *  @param length length of the nucleotide sequence
 *  @return pointer to char 
 * 
 *  Encodes a sequence of plain nucleotides
 */
char* encodeBases(char *dest, char* src, unsigned int length);

/**
 *  @brief Decodes a sequence of encoded nucleotides
 *  @param dest pointer to destination char with plain nucleotides
 *  @param src pointer to char with encoded nucleotides
 *  @param length length of the nucleotide sequence
 *  @return pointer to char 
 * 
 *  Decodes a sequence of encoded nucleotides
 */
char* decodeBases(char *dest, char* src, unsigned int length);

extern int nA;
extern int AA, CC, GG, TT;

extern int table[128];
extern int rev_table[4];

void initReplaceTable_bs(const char *str);

#endif	/*    STRING_UTILS_H	*/
