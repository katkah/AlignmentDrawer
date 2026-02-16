#include<iostream>
#include<fstream>
#include<vector> 
#include<iterator>
#include<sstream>
#include<cassert>
#include<Board.h>
#include<array>
#include<math.h>
using namespace std;
using namespace LibBoard;

const int free_space = 10;
const int font_size_1 = 60;
const int font_size_2 = 90;
const int font_size_3 = 100;
const int font_size_4 = 8;
const int free_space_2 = font_size_3 * 2;
const int len_rect1 = 10;
const double len_char = 60.1545; //len_char for font_size_3 = 100
const double max_x = 900 * len_char;
Color DarkBlue = Color( 58, 95, 205 ); 
Color DarkRed = Color( 238, 0, 0 ); 
Color DarkGreen = Color( 0, 100, 0 ); 
Color DarkGray = Color( 64, 64, 64 );

double img_len( string seq ) {
    return seq.length() * len_char;
}

int number_of_equal_elements(string::iterator seq_start, string::iterator  seq_end, vector<float>::iterator entropy_start){
    int num_of_equal = 1;
    while((seq_start + 1) < seq_end  && *(seq_start + 1) != '-'){
        int color1 = 256 - *(entropy_start) * 100;
        int color2 = 256 - *(entropy_start + 1) * 100;
        if (color1 == color2){
            ++num_of_equal; 
            ++seq_start;
            ++entropy_start;
        }else return num_of_equal;
    }
    
    return num_of_equal;
}

vector<pair<int,int>> print_deletion(const string &seq){
    vector<pair<int,int>> deletion_pos;
    unsigned i = 0;
    while (i < seq.size()) {
        if(seq[i] == '-'){
            unsigned start = i;
            while (i < seq.size() && seq[i] == '-')
                ++i;
            deletion_pos.emplace_back(start, i);
        }
        ++i;
    }
    return deletion_pos;
}

void print_letter(string letter, int X_pos, int Y_pos, Board &board){
    switch (letter[0]){
        case 'A':
            board << Text( X_pos, Y_pos, "A", Fonts::Font( Fonts::Monospace ), font_size_3, Color::Green );
            break;
        case 'C':
            board << Text( X_pos, Y_pos, "C", Fonts::Font( Fonts::Monospace ), font_size_3, DarkBlue );
            break;
        case 'G':
            board << Text( X_pos, Y_pos, "G", Fonts::Font( Fonts::Monospace ), font_size_3, Color::Black );
            break;
        case 'T':
            board << Text( X_pos, Y_pos, "T", Fonts::Font( Fonts::Monospace ), font_size_3, DarkRed );
            break;
        default:
            board << Text( X_pos, Y_pos, letter, Fonts::Font( Fonts::Monospace ), font_size_3, Color::Black );
            break;
    }
}

void print_heads(int X_pos, int Y_pos, Board &board, const vector<string> heads, const vector<int> heads_len){
    for( int j = 0; j < heads.size(); j++ ){ //writes heads
        int Y = Y_pos + j * (-(font_size_3 + free_space));
        board << Text( 0, Y, heads[j], Fonts::Font(Fonts::Monospace ), font_size_3, DarkGray );
    }
}
string multiply( string unit, int times ) {
    string ret;
    while ( times-- )
        ret += unit;
    return ret;
}


int max_element( vector <int> nums) {//finds max_element in vector
    
    int max = nums[0];
    for( int i = 1; i < nums.size(); i++ ){
        if (max < nums[i]) max = nums[i];
    }
    return max;
}


void print_coordinate(stringstream &seq_coordinate, int X, int Y, Board &board){
    Y = Y - ( font_size_3 + free_space );
    board << Line( X + len_char/2, Y + font_size_1, 
    X + len_char/2, Y + 1.5 * font_size_1, 
    DarkGray, len_char/8 );//marks for sequence coordinates 
    board << Text( X, Y, seq_coordinate.str(), Fonts::Font( Fonts::Monospace ), font_size_1, DarkGray );
}

vector<int> compute_borders( vector<string> sequences ) {
    vector<int> borders;
    
    for(string seq : sequences) {
        if (seq.size() != sequences[0].size()) {
            cerr << "sequences are not the same length"  << endl;
            exit( 1 );    
        }    
    }
    bool prevSame = true;
    for(int i = 0; i < sequences[0].size(); i++){
        int thisSame = true;
        for(int j = 0; j < sequences.size(); j++){
            if(sequences[0][i] != sequences[j][i]){
                thisSame = false;
            }
        }
        if(thisSame != prevSame){
            borders.push_back(i);
            prevSame = !prevSame;
        }
    }
    borders.push_back( sequences[0].size() );
    return borders;
}

int main( int argc, char** argv ){
    Board board;

    string line;
        
    vector<string> sequences;
    vector<string> heads;
    vector<int> borders;

    while (cin.good()){
        getline (cin, line);
        if (line.size() == 0){
            continue;
        }else if (line[0] == '>' ) {
            heads.push_back(line);
            sequences.push_back("");    
        } else
            sequences.back()+=line;
    }
    if (sequences.empty()){return 1;}

    borders = compute_borders( sequences );
    
    vector<int> heads_len;
    for( int j = 0; j < sequences.size(); j++ ){ //counts length of heads
        heads_len.push_back( img_len( heads[j] ) );    
    }
    
    


    const int lowest_Y = ( -( font_size_3 + free_space ) ) * sequences.size();
    int break_count;
    int X_pos;
    int Y_pos;
    for( int j = 0; j < sequences.size(); j++ ){//printes sequences
        break_count = 1;
        X_pos = 0 + max_element( heads_len ) + len_char;
        Y_pos = j * ( -( font_size_3 + free_space ) );
        bool equal = true;
        if ( borders[0] == 0 ) equal = false;
        for ( int i = (equal ? 0 : 1); i < borders.size(); i++ ) {
            if ( !equal ) {//printes block of letters
                string sub = sequences[j].substr( borders[i-1], borders[i]-borders[i-1] );
                stringstream seq_coordinate;
                seq_coordinate << ( borders[i-1] + 1 ); //absolute sequence coordinates
                
                for( int k = 0; k < sub.size(); k++ ){
                    if( X_pos >= max_x ){ //breaks lines
                        X_pos = max_element(heads_len) + len_char;
                        Y_pos += (lowest_Y - free_space_2);
                        break_count ++;
                    }
                    print_letter(sub.substr(k,1), X_pos,Y_pos,board);
                    if((k == 0) && (j == sequences.size()-1)){                   
                        print_coordinate(seq_coordinate, X_pos, Y_pos, board);
                    }
                    X_pos += len_char;
                }
            } else {//printes boxes
                stringstream num_box;

                if (i == 0) {
                    num_box << borders[0];
                } else {
                    num_box << borders[i]-borders[i-1];
                }
                int num_len = img_len( num_box.str() );
                int len_rect = num_len + len_char;
                if( (X_pos + len_rect) >= max_x){//breaks lines
                    X_pos = max_element(heads_len) + len_char;
                    Y_pos += (lowest_Y - free_space_2);
                    break_count++;
                }
                Rectangle b( X_pos, Y_pos + font_size_3, len_rect, font_size_3, Color::Green, DarkGreen, 5 );
                board << b;
                board << Text( X_pos + ( len_rect - num_len ) / 2, Y_pos + 5, num_box.str(),
                 Fonts::Font( Fonts::Monospace ), font_size_2, Color::White );
                X_pos += len_rect;
            }
            equal = !equal;
                
        }
    }
    
    board << Line( max_x + 3 * len_char, Y_pos, max_x + 3 * len_char, Y_pos, Color::Black, 22 );//marks for sequence coordinates
    int Y = 0;
    for( int i = 1; i <= break_count; i++){//prints heads
        print_heads( 0, Y, board, heads, heads_len);
        Y += (lowest_Y - free_space_2);
    }
    
    
    //board.saveSVG( NULL );

    /******************* second image
    ****************************************************************************************/
    
    Board board1;
    using fourtuple = array<int, 4>;
    vector<fourtuple> occurances;
    vector<float> entropy;
    int a,c,g,t,delece;
    for(int i = 0; i < sequences[0].size(); i++){
        a = c = g = t = delece = 0;
        for(int j = 0; j < sequences.size(); j++){
            switch ( sequences[j][i] ) {
                case 'A' :
                    a++;
                    break;
                case 'C' :
                    c++;
                    break;
                case 'G' :
                    g++;
                    break;
                case 'T' :
                    t++;
                    break;
                case '-' :
                    delece++;
                    break;
            }
        }
        int sum = sequences.size() - delece;
        fourtuple occ = {{a,c,g,t}};
        occurances.push_back(occ);
        entropy.push_back( 0 );
        for ( int k = 0; k < 4; ++k) {//counts entropy of a column
            if(occ[k] != 0) {
                float freq = float(occ[k]) / sum;
                entropy.back() -= freq * log2( freq );

            }
        }
    }
    Y = 0;
    print_heads( 0, Y, board1, heads, heads_len);
    int color = 0;
    int n_equal_elements = 1;
    for( int s = 0; s < sequences.size(); ++s ){//printes sequences
        X_pos = 0 + max_element( heads_len ) + len_char;
        Y_pos = s * ( -( font_size_3 + free_space ) );
        int position = 0;
        auto deletion_pos = print_deletion(sequences[s]);
        auto deletion_iterator = deletion_pos.begin();
        auto seq_start_iterator = sequences[s].begin();
        auto seq_end_iterator = sequences[s].end();
        auto entropy_start_iterator = entropy.begin();
        while (position < sequences[0].size()){
            if ( sequences[s][position] != '-' ){
                color =  256 - entropy[position] * 100;//100 natahuje barevné spektrum 
                auto rect_color = Color::Color(color, color, color);
                n_equal_elements = number_of_equal_elements(seq_start_iterator, seq_end_iterator, entropy_start_iterator);
                Rectangle r( X_pos, Y_pos + font_size_3, len_rect1 * n_equal_elements, font_size_3, rect_color, rect_color, 0 );
                board1 << r;
                X_pos += len_rect1 * n_equal_elements;
                seq_start_iterator += n_equal_elements;
                entropy_start_iterator += n_equal_elements;
                position += n_equal_elements;
            }else{
                int deletion_length = len_rect1 * (deletion_iterator->second - deletion_iterator->first);
                Rectangle r( X_pos, Y_pos + font_size_3/2, deletion_length, font_size_4, Color::Black , Color::Black, 0 );
                board1 << r;
                X_pos += deletion_length;
                position = deletion_iterator->second;
                ++deletion_iterator;
            }
        }
    }
    board1.saveSVG( NULL ); 
//weblogo pro uživatelskou parametrizaci
    return 0;
}
