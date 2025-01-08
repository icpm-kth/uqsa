#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <Rdefines.h>
typedef SEXP Rdata;

struct stack {
	char *buffer;
	size_t size;
	size_t n;
};

struct stack* stack_alloc(size_t N){
	struct stack *s = malloc(sizeof(struct stack));
	s->size = N;
	s->n = 0;
	s->buffer = malloc(sizeof(char)*s->size);
	return s;
}

void stack_free(struct stack **s){
	free((*s)->buffer);
	(*s)->n = 0;
	(*s)->size = 0;
	*s=NULL;
}

void push(struct stack *s, char val){
	if (s==NULL) return;
	if (s->n == s->size - 1){
		s->size += 12;
		s->buffer = realloc(s->buffer,s->size);
	}
	s->buffer[s->n++] = val;
	s->buffer[s->n] = '\0'; /* for safe printing */
}

char pop(struct stack *s){
	char c='\0';
	if (s->n > 0) {
		c = s->buffer[--(s->n)];
		s->buffer[s->n] = '\0';
	}
	return c;
}

int empty(struct stack *s){
	return s->n == 0;
}

void stack_print(struct stack *s){
	printf("«%s»\n",s->buffer);
}

Rdata replace_pow(Rdata str){
	const char *input = CHAR(STRING_ELT(str,0));
	char *ptr = input;
	char *hat; /* points to U+005E ^ CIRCUMFLEX ACCENT if found in str */
	char *open=NULL, *close=NULL; /* first opening parenthesis related to ^ */
	struct stack *s = stack_alloc(24);
	Rdata w = PROTECT(NEW_CHARACTER(1));
	char *out;
	int n_prefix,n_base,n_exponent,n_suffix;
	hat=strchr(ptr,'^');
	if (hat) {
		ptr = hat-1;
		while (ptr > input && open==NULL){
			switch (*ptr){
			case ')': push(s,')'); break;
			case '(': pop(s); break;
			case '-':
			case '+':
			case '*':
			case '/':
			case ' ':
				if (empty(s)){
					open = ptr+1;
				}
			}
			ptr--;
		}
		n_base = hat - open;
		n_exponent = close - hat;
		n_prefix = open - input;
		n_suffix = strlen(close);
		out = malloc(n_prefix + n_base + n_exponent + n_suffix + 4);
		memcpy(out,input,n_prefix);
	}
	SET_STRING_ELT(w, 0, mkChar("a custom string constant"));
	stack_free(&s);
	UNPROTECT(1);
	return w;
}
