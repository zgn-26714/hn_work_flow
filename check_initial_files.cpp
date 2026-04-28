#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

void change_inp(std::string inp_file, std::string pdb_name, std::string idx_str);
void change_pdb_name(std::string pdb_file, std::string name);
std::string find_pdb_file(const std::vector<std::string>& pdb_files, const std::string& mol_name);

// 按空格分割字符串
std::string trim(const std::string& s)
{
    size_t start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

std::vector<std::string> split_by_space(const std::string& str)
{
    std::vector<std::string> result;
    std::istringstream iss(str);
    std::string token;
    while (iss >> token)
    {
        result.push_back(token);
    }
    return result;
}

// 判断top文件中是否包含目标字符串
std::string find_in_top(const std::string& filename, const std::string& target)
{
    std::ifstream fin(filename);
    if (!fin.is_open())
    {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return "";
    }

    std::string line;
    bool in_molecules = false;

    while (std::getline(fin, line))
    {
        line = trim(line);
        if (line.empty()) continue;
        if (line[0] == ';') continue;
        
        if (!in_molecules)
        {
            if (line == "[ molecules ]")
                in_molecules = true;
            continue;
        }

        // 如果已经在 [ molecules ] 段中，遇到下一个 section 就停止
        if (line.front() == '[' && line.back() == ']')
            break;

        // 去掉行尾注释
        size_t comment_pos = line.find(';');
        if (comment_pos != std::string::npos)
            line = trim(line.substr(0, comment_pos));

        if (line.empty())
            continue;

        // 读取第一列：分子名
        std::istringstream iss(line);
        std::string mol_name;
        iss >> mol_name;

        if (mol_name == target)
            return line;
    }

    return "";
}

void get_file_name(std::string &top_file, std::string &inp_file,
                   std::vector<std::string>& pdb_file,
                   std::string TOPorBulkTop){
    const char* top_env, *packmol_env;
    const char* top_var_name, *packmol_var_name;
    if(TOPorBulkTop == "TOP"){
        top_env = std::getenv("TOP");
        packmol_env = std::getenv("packmol");
        top_var_name = "TOP";
        packmol_var_name = "packmol";
    }
    else{
        top_env = std::getenv("bulk_top");
        packmol_env = std::getenv("bulk_pa");
        top_var_name = "bulk_top";
        packmol_var_name = "bulk_pa";
    }
    if (top_env != nullptr){
        top_file = std::string(top_env) + ".top";
    }
    else{
        std::cerr << "Environment variable '" << top_var_name << "' is not set." << std::endl;
        std::exit(1);
    }

    if (packmol_env != nullptr){
        inp_file = std::string(packmol_env) + ".inp";
    }
    else{
        std::cerr << "Environment variable '" << packmol_var_name << "' is not set." << std::endl;
        std::exit(1);
    }

    std::ifstream fin(top_file);
    if (!fin.is_open())
    {
        std::cerr << "Error: cannot open top file: " << top_file << std::endl;
        std::exit(1);
    }
    fin.close();

    /*get pdb file*/
    pdb_file.clear();
    fin.open(inp_file);

    if (!fin.is_open())
    {
        std::cerr << "Error: cannot open input file " << inp_file << std::endl;
        std::exit(1);
    }
    std::string line;
    while (std::getline(fin, line))
    {
        std::istringstream iss(line);
        std::string key, filename;
        iss >> key >> filename;

        // 只处理以 "structure" 开头的行
        if (key == "structure" && !filename.empty())
        {
            pdb_file.push_back(filename);
        }
    }
    fin.close();
}


int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " \"name1 name2 name3\"" << std::endl;
        return 1;
    }

    const char* mol_env = std::getenv("MOL_name");
    if (mol_env == nullptr)
    {
        std::cerr << "Error: Environment variable 'MOL_name' is not set." << std::endl;
        std::cerr << "Please define MOL_name in the INPUT file." << std::endl;
        return 1;
    }
    std::string MOL_name(mol_env);
    std::vector<std::string> sel_name = split_by_space(MOL_name);

    std::string top_file, inp_file;
    std::vector<std::string> pdb_files;

    int begin, end;
    std::string mode = argv[1];

    if (mode == "all")
    {
        begin = 0;
        end = 2;
    }
    else if (mode == "top")
    {
        begin = 0;
        end = 1;
    }
    else if (mode == "bulk")
    {
        begin = 1;
        end = 2;
    }
    else
    {
        std::cerr << "Error: invalid mode '" << mode
                  << "'. Use all, top, or bulk." << std::endl;
        return 1;
    }


    std::vector<std::string> set_rules = {"TOP", "bulk"};

    

    for (const auto& name : sel_name)
    {
        std::string pdb_name;
        for (int i = begin; i < end; ++i){
            get_file_name(top_file, inp_file, pdb_files, set_rules[i]);
            
            // std::cout << "===== File info for set_rule: " << set_rules[i] << " =====" << std::endl;
            // std::cout << "Top file: " << top_file << std::endl;
            // std::cout << "Input file: " << inp_file << std::endl;

            // std::cout << "ITP files: ";
            // for (const auto& itp : itp_files)
            //     std::cout << itp << " ";
            // std::cout << std::endl;

            // std::cout << "PDB files: ";
            // for (const auto& pdb : pdb_files)
            //     std::cout << pdb << " ";
            // std::cout << std::endl;

            // std::cout << "=========================================" << std::endl;
            // std::cout<<"find mol: "<<name<<std::endl;
            std::string idx_top = find_in_top(top_file, name);//mol + num
            if (!idx_top.empty())
            {   
                //另一种思路是根据top文件中的mol num找到inp文件中对应的structure行，找到对应的pdb文件，再根据pdb文件中的分子名找到对应的pdb文件
                //但是可行性需要验证（当top中的名字和pdb中的，itp中的名字不一致会发生什么？）
                pdb_name = find_pdb_file(pdb_files, name);//找到名字对应的pdb文件，没找到会报错退出
                //std::cout<<"get_pdb_name: "<<pdb_name<<std::endl;
                change_inp(inp_file, pdb_name, idx_top);//改变inp文件中对应的pdb文件的名字,并把下一行的number改成和top中一致

                
                // break; // 找到了就查下一个 name
            }
            else{
                std::cerr << "Error: molecule name '" << name << "' not found in top file." << std::endl;
                std::cerr << "Please check your INPUT or  check your initial file by yourself!" << std::endl;
                std::exit(1);
            }

            // std::string idx_inp = find_in_inp(inp_file, name);
            // if (idx_inp.c_str() != nullptr)
            // {   
            //     found = true;
            //     change_top<std::string>(top_file, idx_inp);
            //     change_itp(itp_file, idx_inp);
            //     change_pdb(pdb_file, idx_inp);
            //     break;  // 找到了就查下一个 name
            // }


    //itp中的分子名必然和pdb中的一致，所以不需要单独查itp文件了
            // std::string idx_itp = find_in_itp(itp_file, name);
            // if (idx_itp.c_str() != nullptr)
            // {   
            //     found = true;
            //     change_top(top_file, idx_itp);
            //     change_inp(inp_file, idx_itp);
            //     change_pdb(pdb_file, idx_itp);
            //     break;  // 找到了就查下一个 name
            // }

            // std::string idx_pdb = find_in_pdb(pdb_file, name);
            // if (idx_pdb.c_str() != nullptr)
            // {   
            //     found = true;
            //     change_top(top_file, idx_pdb);
            //     change_inp(inp_file, idx_pdb);
            //     change_itp(itp_file, idx_pdb);
            //     break;  // 找到了就查下一个 name
            // }
            // if(!found){
            //     std::cerr << "Error: cannot find molecule name '" << name << "' in any of the files." << std::endl;
            //     std::exit(1);
            // }
        }
        change_pdb_name(pdb_name, name);//重命名pdb文件
    }
    

    return 0;
}

void change_inp(std::string inp_file, std::string pdb_name, std::string idx_str){
    std::vector<std::string> ref_str = split_by_space(idx_str);
    std::string ref_mol_name = ref_str[0];
    int ref_mol_num = std::stoi(ref_str[1]);  
    std::ifstream fin(inp_file);
    if (!fin.is_open())
    {
        std::cerr << "Cannot open file: " << inp_file << std::endl;
        std::exit(1);
    }

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(fin, line))
    {
        lines.push_back(line);
    }
    fin.close();

    int structure_idx = -1;
    int match_count = 0;

    for (int i = 0; i < static_cast<int>(lines.size()); ++i)
    {
        std::string cur = trim(lines[i]);
        std::istringstream iss(cur);

        std::string key, value;
        iss >> key >> value;

        if (key == "structure" && value == pdb_name)
        {
            ++match_count;
            structure_idx = i;
        }
    }

    if (match_count == 0)
    {
        std::cerr << "Error: structure " << pdb_name
                << " not found in " << inp_file << std::endl;
        std::exit(1);
    }

        int number_idx = -1;

    // 在对应的 structure 块里继续向下找 number
    for (int j = structure_idx + 1; j < static_cast<int>(lines.size()); ++j)
    {
        std::string sub = trim(lines[j]);
        if (sub.empty()) continue;

        std::istringstream jss(sub);
        std::string subkey;
        jss >> subkey;

        if (subkey == "end")
        {
            std::string sub2;
            jss >> sub2;
            if (sub2 == "structure")
                break;
        }

        if (subkey == "number")
        {
            int old_num;
            if (!(jss >> old_num))
            {
                std::cerr << "Error: failed to parse number in structure block of "
                          << pdb_name << std::endl;
                std::exit(1);
            }
            number_idx = j;
            break;
        }
    }

    if (number_idx == -1)
    {
        std::cerr << "Error: cannot find 'number' in structure block of "
                  << pdb_name << " in " << inp_file << std::endl;
        std::exit(1);
    }
    // 修改 structure 那一行
    lines[structure_idx] = "structure " + ref_mol_name + ".pdb";
    lines[number_idx] = "number " + std::to_string(ref_mol_num);
    std::ofstream fout(inp_file);
    if (!fout.is_open())
    {
        std::cerr << "Cannot write file: " << inp_file << std::endl;
        std::exit(1);
    }
    for (const auto& l : lines)
    {
        fout << l << "\n";
    }
    fout.close();
   
}

std::string find_pdb_file(const std::vector<std::string>& pdb_files, const std::string& mol_name)
{
    for (const auto& pdb_file : pdb_files)
    {
        std::ifstream fin(pdb_file);
        if (!fin.is_open())
        {
            std::cerr << "Cannot open file: " << pdb_file << std::endl;
            continue;
        }

        std::string line;
        while (std::getline(fin, line))
        {
            if (line.compare(0, 4, "ATOM") != 0) continue;

            std::vector<std::string> fields = split_by_space(line);
            if (fields.size() >= 4 && fields[3] == mol_name)
            {
                return pdb_file;
            }
        }
    }

    std::cerr << "Error: molecule name '" << mol_name
              << "' not found as the 4th field of any ATOM record in ";
    for (const auto& pdb_file : pdb_files){
        std::cerr << pdb_file << " ";
    }
    std::cerr << std::endl;

    std::exit(1);
}

void change_pdb_name(std::string pdb_file, std::string name){
    std::string new_file = name + ".pdb";
    std::string cmd = "mv \"" + pdb_file + "\" \"" + new_file + "\"";

    if(pdb_file == new_file){
        return;
    }
    int ret = std::system(cmd.c_str());
    if (ret != 0)
    {
        std::cerr << "Error: failed to rename file from "
                  << pdb_file << " to " << new_file << std::endl;
        std::exit(1);
    }
}