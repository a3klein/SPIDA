const yearSpan = document.querySelector("[data-year]");
const docLinks = document.querySelector("#doc-links");
const docContent = document.querySelector("#doc-content");

if (yearSpan) {
	yearSpan.textContent = new Date().getFullYear();
}

const docs = [
	{ label: "README", file: "README.md" },
	{ label: "Quick start", file: "quick_start.md" },
	{ label: "Configuration", file: "configuration.md" },
	{ label: "Software setup", file: "software_setup.md" },
	{ label: "Deconvolution", file: "deconvolution.md" },
];

const getDocFromHash = () => {
	const raw = window.location.hash.replace("#", "");
	if (!raw) {
		return docs[0].file;
	}

	const match = docs.find((doc) => doc.file.toLowerCase() === raw.toLowerCase());
	return match ? match.file : docs[0].file;
};

const renderLinks = () => {
	if (!docLinks) {
		return;
	}

	docLinks.innerHTML = "";
	for (const doc of docs) {
		const item = document.createElement("li");
		const link = document.createElement("a");
		link.href = `#${doc.file}`;
		link.textContent = doc.label;
		link.classList.add("doc-link");
		item.appendChild(link);
		docLinks.appendChild(item);
	}
};

const highlightActiveLink = (file) => {
	const links = document.querySelectorAll(".doc-link");
	for (const link of links) {
		link.classList.toggle("active", link.getAttribute("href") === `#${file}`);
	}
};

const renderMarkdown = (markdown) => {
	if (!docContent) {
		return;
	}

	if (window.marked) {
		docContent.innerHTML = window.marked.parse(markdown, { mangle: false });
	} else {
		docContent.textContent = markdown;
	}
};

const loadDoc = async (file) => {
	try {
		const response = await fetch(file);
		if (!response.ok) {
			throw new Error(`Failed to load ${file}`);
		}
		const markdown = await response.text();
		renderMarkdown(markdown);
		highlightActiveLink(file);
	} catch (error) {
		renderMarkdown(`## Unable to load ${file}\n\n${error.message}`);
	}
};

const init = () => {
	renderLinks();
	const file = getDocFromHash();
	loadDoc(file);
};

window.addEventListener("hashchange", () => {
	const file = getDocFromHash();
	loadDoc(file);
});

window.addEventListener("DOMContentLoaded", init);
