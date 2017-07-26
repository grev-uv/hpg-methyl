#include <stdio.h>
#include <string.h>

#include "../log.h"

int main(void) {

	printf("Testing equals(casa, casa)\n");
	printf("Result: %i\n\n", equals("casa", "casa"));

	printf("Testing equals(casa, saca)\n");
	printf("Result: %i\n\n", equals("casa", "saca"));

	printf("Testing equals_ignore_case(casa, casa)\n");
	printf("Result: %i\n\n", equals_ignore_case("casa", "casa"));

	printf("Testing equals_ignore_case(casa, cAsa)\n");
	printf("Result: %i\n\n", equals_ignore_case("casa", "cAsa"));

	printf("Testing equals_ignore_case(_ca/sa, _Ca/Sa)\n");
	printf("Result: %i\n\n", equals_ignore_case("_ca/sa", "_Ca/Sa"));

	printf("Testing is_numeric(1456)\n");
	printf("Result: %i\n\n", is_numeric("1456"));

	printf("Testing is_numeric(45t6)\n");
	printf("Result: %i\n\n", is_numeric("45t6"));



	printf("Testing to_lower_case(Lap_TOp*)\n");
	char *str1 = (char*)malloc(20*sizeof(char));
	strcpy(str1, "Lap_TOp*");
	printf("Result: %s\n\n", to_lower_case(str1));

	printf("Testing to_upper_case(Lap_TOp*)\n");
	char str2[] = "Lap_TOp*";
	printf("Result: %s\n\n", to_upper_case(str2));



	char *str3 = (char*)malloc(20*sizeof(char));
	char *str4 = (char*)malloc(20*sizeof(char));
	printf("Testing chop(Lap_TOp*)\n");
	strcpy(str3, "Lap_TOp*");
	chop(str3);
	printf("Result: %s\n\n", str3);

	printf("Testing chop_at(Lap_TOp*)\n");
	strcpy(str3, "Lap_TOp*");
	int str3_len = strlen(str3);
	chop_at(str3, str3_len-1);
	printf("Result: %s\n\n", str3);

	printf("Testing chomp(Lap_TOp*\\n)\n");
	strcpy(str3, "Lap_TOp*\n");
	chomp(str3);
	printf("Result: %s\n\n", str3);

	printf("Testing chomp_at(Lap_TOp*\\n)\n");
	strcpy(str3, "Lap_TOp*\n");
	str3_len = strlen(str3);
	chomp_at(str3, str3_len-1);
	printf("Result: %s\n\n", str3);



	printf("Testing remove_char(Lap_TOpa*, 'a')\n");
	strcpy(str3, "Lap_TOpa*");
	remove_char(str3, 'a');
	printf("Result: %s\n\n", str3);

	printf("Testing remove_char_at(Lap_TOpa*, 4)\n");
	strcpy(str3, "Lap_TOpa*");
	remove_char_at(str3, 4);
	printf("Result: %s\n\n", str3);

	printf("Testing remove_str(Laap_TOpaa*aa, aa)\n");
	strcpy(str3, "Laap_TOpaa*aa");
	strcpy(str4, "aa");
	remove_str(str3, str4);
	printf("Result: %s\n\n", str3);

	printf("Testing remove_start(Lap_TOpa*, 3)\n");
	strcpy(str3, "Lap_TOpa*");
	remove_start(str3, 3);
	printf("Result: %s\n\n", str3);

	printf("Testing remove_end(Lap_TOpa*, 4)\n");
	strcpy(str3, "Lap_TOpa*");
	remove_end(str3, 4);
	printf("Result: %s\n\n", str3);



	printf("Testing trim( \t Lap_TOpa* \\t )\n");
	strcpy(str3, " \t Lap_TOpa* \t ");
	trim(str3);
	printf("Result: '%s'\n\n", str3);

	printf("Testing ltrim2( \\n \\t Lap_TOpa* \\t)\n");
	strcpy(str3, " \n \t Lap_TOpa* \t");
	ltrim2(str3);
	printf("Result: '%s'\n\n", str3);

	printf("Testing rtrim2(Lap_TOpa* \\t )\n");
	strcpy(str3, " Lap_TOpa* \t ");
	rtrim2(str3);
	printf("Result: '%s'\n\n", str3);

	printf("Testing strip( \t Lap_TOpa* \\t )\n");
	strcpy(str3, " \t Lap_TOpa* \t ");
	strip(str3);
	printf("Result: '%s'\n\n", str3);

	printf("Testing lstrip(  \\t Lap_TOpa* \\t)\n");
	strcpy(str3, "  \t Lap_TOpa* \t");
	lstrip(str3);
	printf("Result: '%s'\n\n", str3);

	printf("Testing rstrip( Lap_TOpa* \\t )\n");
	strcpy(str3, " Lap_TOpa* \t ");
	rstrip(str3);
	printf("Result: '%s'\n\n", str3);

	return 1;
}
