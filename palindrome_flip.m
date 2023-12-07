function palindrome = palindrome_flip(sequence)
s = size(sequence);
flipped_seq = fliplr(sequence);
palindrome = flipped_seq;
for i = 1:s(2)
    if flipped_seq(i)=='A'
        palindrome(i) = 'T';
    elseif flipped_seq(i)=='T'
        palindrome(i) = 'A';
    elseif flipped_seq(i)=='G'
        palindrome(i) = 'C';
    elseif flipped_seq(i)=='C'
        palindrome(i) = 'G';
    else
        disp('Invalid base in file.')
    end
end