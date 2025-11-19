def hamming_distance(str1, str2):
    """
    Calculate the Hamming distance between two strings.
    If the strings are not the same length, pad the shorter one with underscores ('_').
    """
    # Make both strings the same length by padding
    max_len = max(len(str1), len(str2))
    str1 = str1.ljust(max_len, "_")
    str2 = str2.ljust(max_len, "_")

    # Compute the number of differing characters
    distance = sum(ch1 != ch2 for ch1, ch2 in zip(str1, str2))
    return distance


# Example usage
slack_username = "CHRISTOPHER"
x_handle = "CHRISSEUN11"

distance = hamming_distance(slack_username, x_handle)
print(f"Hamming distance between '{slack_username}' and '{x_handle}' is: {distance}")
