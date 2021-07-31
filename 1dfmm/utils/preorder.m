function pre_ord = preorder(i, tr, ch)

% return the pre order traversal of subtree rooted at i

if isempty(ch{i})
    pre_ord = i;
    return
else
    c1 = ch{i}(1);
    c2 = ch{i}(2);
    pre_ord_c1 = preorder(c1, tr, ch);
    pre_ord_c2 = preorder(c2, tr, ch);
    pre_ord = [i, pre_ord_c1, pre_ord_c2];
    return
end