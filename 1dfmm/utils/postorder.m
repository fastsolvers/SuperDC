function post_ord = postorder(i, tr, ch)

% return the post order traversal of subtree rooted at i

if isempty(ch{i})
    post_ord = i;
    return
else
    c1 = ch{i}(1);
    c2 = ch{i}(2);
    post_ord_c1 = postorder(c1, tr, ch);
    post_ord_c2 = postorder(c2, tr, ch);
    post_ord = [post_ord_c1, post_ord_c2, i];
    return
end