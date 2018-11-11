font.match {
    color: red;
}

div.gather {
    display: none;
}

div.gatherdetail {
    display: none;
}

input[type=checkbox].gather:checked ~ div.gather {
    display: inline-block;
}

input[type=checkbox].gatherdetail:checked ~ div.gatherdetail {
    display: inline-block;
}

