void test(){

  string str1 = "aho";
  string str2 = "tyou aho";
  string str3 = "baka";

  if(strstr(str2.c_str(),str1.c_str()) != NULL) cout << str1 << " is included in " << str2 << endl;
  if(strstr(str3.c_str(),str1.c_str()) != NULL) cout << str1 << " is included in " << str3 << endl;;
  if(strstr(str2.c_str(),str1.c_str()) == NULL) cout << str1 << " is not included in " << str2 << endl;
  if(strstr(str3.c_str(),str1.c_str()) == NULL) cout << str1 << " is not included in " << str3 << endl;;

}
