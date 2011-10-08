EncNGram::EncNGram(unsigned char* dataref, char size) {
       //create a copy of the character data (will take less space than storing pointers anyhow!)
       _size = size;
       data = new unsigned char[size];
       n = 1;
       for (int i = 0; i < size; i++) {
            data[i] = dataref[i];
       }
}

const char EncNGram::size() {
    return _size;
}

const int EncNGram::n() {
    int count = 1; 
    for (int i = 0; i < _size; i++) {
        if (data[i] == 0) count++;
    }    
    return count;
}

EncNGram::EncNGram slice(const int begin,const int length) {
    return getencngram(begin, length, data, _size);
}


EncNGram getencngram(const int index, const int n, unsigned char *line, const int size) {
    int currentindex = 0;
    int beginpos = -1;
    int endpos = -1;
    for (int i = 0; i < size; i++) {
        if (line[i] == 0) {
            currentindex++;
            if (currentindex == index) {
                beginpos = i+1;
            } else if (currentindex == index + n) {
                endpos = i - 1;
            }
        }        
    }
    if (endpos == -1) {
        endpos = size - 1;
    }
    const char size = (char) (endpos - beginpos + 1);
    return EncNGram(line + pos, size); 
}
